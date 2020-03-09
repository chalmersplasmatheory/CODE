classdef SmartLUSolver < Solver

    properties (SetAccess = protected)
        timeIndex = 1;
        fattimeIndex = [];
        nSaves = 30;
        saveIndices;
        factor_L;
        factor_U;
        factor_P;
        factor_Q;
    end
    
    properties (Constant)
        VERSION = 1.0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%   Constructor    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function this = SmartLUSolver(state,eqSettings)
            this@Solver(state,eqSettings)
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% equation solving  %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        %if output is supplied, this will be modified, usefull when
        %extending the distributions
        %then remember that this is pointers so you dont need to collect
        %the output
        function output =  takeTimeSteps(this,output)
            if nargin < 2
                output = Output();
            end
            this.initsaveIndices()
            fbefore = [];
            for i = this.timeIndex:(this.state.timeGrid.nTimeSteps-1)
                rhs = this.fattimeIndex;
                if i > 1 %only to compare with CODE
                    rhs(1:this.state.momentumGrid.Ny) = rhs(1:this.state.momentumGrid.Ny) + this.enforceParticleAndHeat;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%% Explicit Operators %%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for EO = this.explicitOperators
                    EO{1}.generateOperatorMatrix(i);
                    rhs = rhs + this.state.timeGrid.dts(i)*EO{1}.getMatrix*this.fattimeIndex;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%% Sources %%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for S = this.sources
                    rhs = rhs + this.state.timeGrid.dts(i)*S{1}.getSourceVec(this.fattimeIndex,i);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%% Implicit Operators %%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                impMat = sparse(this.state.momentumGrid.matrixSize, this.state.momentumGrid.matrixSize);
                matrixChange = zeros(size(this.implicitOperators));
                k = 1;
                for IO = this.implicitOperators
                    %matrixChange(k) = IO{1}.generateOperatorMatrix(i+1); %this is correct but code does below
                    matrixChange(k) = IO{1}.generateOperatorMatrix(i);
                    k = k+1;
                end
                if any(matrixChange) || this.state.timeGrid.dtsHasChanged(i)
                    for IO = this.implicitOperators
                        if isa(IO{1},'DiagonalPlusBoundaryCondition')
                            impMat = impMat + IO{1}.getMatrix();
                        else
                            impMat = impMat + this.state.timeGrid.dts(i)*IO{1}.getMatrix();
                        end
                    end
                    [this.factor_L, this.factor_U, this.factor_P, this.factor_Q] = lu(impMat);
                    disp('inversion done')
                end
                
                rhs(1) = 0;
                % Handle boundary condition at y = yMax:
                if this.state.momentumGrid.yMaxBoundaryCondition ~= 3
                    rhs((1:this.state.momentumGrid.Nxi)*this.state.momentumGrid.Ny) = 0;
                end
                
                this.timeIndex = i+1;
                fbefore = this.fattimeIndex;
                this.fattimeIndex = this.factor_Q * (this.factor_U \ (this.factor_L \ (this.factor_P * rhs)));
                
                if i == 1 %save first step at first time point, but want to save last step at last time point. To not invert one uneccessary time, make first save special
                    disp('saving')
                    output.save(i,this.state, fbefore,fbefore);
                end
                if any(this.timeIndex == this.saveIndices)
                    disp('saving')
                    output.save(this.timeIndex,this.state, this.fattimeIndex,fbefore);
                end
                
            end
            output.save(this.timeIndex,this.state,this.fattimeIndex,fbefore);
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%     setters        %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %function setf(this,f)
        %    this.fattimeIndex = f;
        %end

        function setf(this,distribution,momentumGrid)
            if ~isa(distribution,'Distribution')
                error('First argument must be a Distribution')
            end
            if nargin == 2
                this.fattimeIndex = distribution.f;
            end
            if nargin == 3
                if ~isa(momentumGrid, 'MomentumGrid')
                    error('Second argument needs to be empty or of type MomentumGrid')
                end
                this.fattimeIndex = interpolateDistribution(distribution,momentumGrid);
            end
        end
        
        % sets so the simulation starts from the time
        function setTime(this,time,timeUnit)
            if nargin < 3
                timeUnit = 'normalized';
            end
            this.timeIndex = find(time>=this.state.timeGrid.gettimesteps(timeUnit),1,'last');
        end
        
        function setnSaves(this,nSaves)
            this.nSaves = nSaves;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%   initializing    %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%   functions       %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        
        function initsaveIndices(this)
            this.saveIndices = unique(round(linspace(...
                                        this.timeIndex,...
                                        this.state.timeGrid.nTimeSteps,...
                                        this.nSaves...
                                    )));
        end
        
        function clearf(this)
            this.fattimeIndex = [];
        end
        
        function setftoZero(this)
            this.fattimeIndex = zeros(this.state.momentumGrid.matrixSize,1);
        end
        
        function resettimeIndex(this)
            this.timeIndex = 1;
        end
        
        function setftoMaxmellian(this)
           maxDis = zeros(this.state.momentumGrid.matrixSize,1);
           maxDis(1:this.state.momentumGrid.Ny) = this.state.physicalParams.nBars(this.timeIndex)...
                                        /this.state.physicalParams.veBars3(this.timeIndex)* ...
                                        exp(...
                                            -this.state.momentumGrid.y2 / ...
                                            this.state.physicalParams.veBars2(this.timeIndex) ...
                                        );
           maxDis(this.state.momentumGrid.Ny) = 0;
           this.fattimeIndex = maxDis;
        end

    end
    
    methods
        %this adds a particle and heat source dependent on what settings are
        %used in eqSetings. These are added to the 0th Legandre Mode.
        
        function rhs0 = enforceParticleAndHeat(this)
            
            Ny = this.state.momentumGrid.Ny;
            y = this.state.momentumGrid.y;
            y2 = this.state.momentumGrid.y2;
            gamma = this.state.momentumGrid.gamma;
            yWeights = this.state.momentumGrid.yWeights;
            energyMoment = y2.*(gamma-1);

            
            f0 = this.fattimeIndex(1:Ny);
            
            nBar = this.state.physicalParams.nBars(this.timeIndex);
            nBarNextStep = this.state.physicalParams.nBars(this.timeIndex+1);
            veBar2NextStep = this.state.physicalParams.veBars2(this.timeIndex+1);
            veBar3 = this.state.physicalParams.veBars3(this.timeIndex);
            veBar2 = this.state.physicalParams.veBars2(this.timeIndex);
            veBar = this.state.physicalParams.veBars(this.timeIndex);
            nInst = (yWeights.*y2) * f0;
            deltaRef2 = this.state.reference.deltaRef^2;
            EOverEc = this.state.physicalParams.EOverEc(this.timeIndex);
            T = this.state.physicalParams.T(this.timeIndex);
            TNextStep = this.state.physicalParams.T(this.timeIndex+1);

            rhs0 = zeros(Ny,1);
             
            %Add a particle source to account for changes in density
            if (abs(this.state.physicalParams.n(this.timeIndex+1)-this.state.physicalParams.n(this.timeIndex))>1e-13 && this.timeIndex ~= this.state.timeGrid.nTimeSteps-1) && ~this.eqSettings.enforceDensityConservation
               %Add a particle source (delta function in time) to account for 
                %the change in density             
                changeFactor = (nBarNextStep-nBar);
                rhs0(1:Ny) = rhs0 + changeFactor /veBar3 ...
                                        * exp(-y2'/veBar2) ;%.* (y2'/veBar2 - 5/2);
            end


            %Enforce heat and/or particle conservation
            if this.timeIndex ~= this.state.timeGrid.nTimeSteps-1 
                %Calculate the source strengths needed
                particleSourceFactor = nBarNextStep - 4/sqrt(pi) * nInst;
                                                             % These operators are not correct for relativistic speeds, so we cut the bulk where the particles have
                                                             % a speed of 0.5 c, when we calculate the temperature. /Albert
                                                             % Ola 2019-04-23: Reconsider energy_mask choice, take as y < 5 instead? 
                                                             % A hot tail seed may now be counted as part of the thermal population!
                energy_mask = y(:)'/veBar < 5;               % Therefore testing this for now.
                fHeat = (yWeights.*energyMoment.*energy_mask) * f0;      
                heatSourceFactor = (deltaRef2*3*sqrt(pi)*nBarNextStep*veBar2NextStep/16 ...
                                            - fHeat)/fHeat; 

                %Add the heat source             
                %   We can add a *dt below if we want. 
                %   If we want to follow the time advance scheme exactly, we
                %   must add a similar term as below for the previous time step 
                %   when using the trapezoid rule.    
                if ~(this.eqSettings.collisionOperator==0 || this.eqSettings.collisionOperator==4)
                    %Turn off the heat source if there is no E-field (useful when
                    %checking conservation properties) and there are no temperature
                    %changes
                    isEField = (EOverEc ~= 0);  
                    if isEField || abs(T-TNextStep) > 1e-13
                        rhs0(1:Ny) = rhs0(1:Ny) + heatSourceFactor*nBar/veBar3 ...
                                        * exp(-y2'/veBar2) .* (y2'/veBar2 - 3/2);
                    end
                end

                %Add the particle source
                if this.eqSettings.enforceDensityConservation
                    %Add a particle source 
                    rhs0(1:Ny) = rhs0(1:Ny) + particleSourceFactor * 1/veBar3 ...
                                        * exp(-y2'/veBar2); % changed particle source to a pure Maxwellian /Ola 2019-04-25
                end

            end
            %rhs0(end) = 0; %is done in CODE I think
        end
        

    end
end

