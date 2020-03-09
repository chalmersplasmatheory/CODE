classdef State <  handle
    
    properties (Constant)
        VERSION = 1.0;
    end
    % %%%%%%%%%%%%%%%%%
    % Properties
    % %%%%%%%%%%%%%%%%%
    %
    % physicalParams - PhysicalParams object shared with timeGrid and reference
    % reference - Reference object shared with all other objects in this class
    % momentumGrid - MomentumGrid object, containing a grid of momentum points
    % timeGrid - TimeGrid object containing timestep vector amongst other

    properties
        physicalParams
        reference
        timeGrid
        momentumGrid
    end
    % %%%%%%%%%%%%%%%%%%%%%
    % These are for the autoInitialGrid
    % %%%%%%%%%%%%%%%%%%%%%
    % NxiScalingFactor - uniformly rescales the predicted Nxi value from 
    %                  - autoInitialGrid by a factor of NxiScalingFactor (default 1)
    % dyBulkScalingFactor - uniformly rescales the desired grid spacing dy
    %                       at y = 0 used by autoInitialGrid (default 1, lower value yields higher resolution)
    % dyTailScalingFactor - uniformly rescales the desired grid spacing dy
    %                       at y = yMax used by autoInitialGrid (default 1, lower value yields higher resolution)
    properties
        NxiScalingFactor    = 1
        Nxi_min             = 3
        Nxi_max             = 150
        pMax_ceiling        = 200
        pMaxIncreaseFactor  = 2
        minPMaxMarginFactor = 4
        dyBulkScalingFactor = 1
        dyTailScalingFactor = 1        
        
        pSwitch = 0
        percentBulk = 0.2
        r = 1.5
    end
    
    methods

        function this = State(TRef,nRef)
            % STATE Initiate a State containing a Reference object shared
            % amongst TimeGrid object, MomentumGrid object and
            % PhysicalParams object. This Reference object ensures that all
            % reference values are consistent between all objects in a
            % State object.
            this.reference = Reference(TRef,nRef);
            this.timeGrid = TimeGrid(this.reference);
            this.momentumGrid = MomentumGrid(this.reference);
            this.physicalParams = PhysicalParams(this.reference,this.timeGrid);
        end
        
        function setInitialRunaway(this, Distribution)
            
        end
        
        function autoInitGrid(this,useScreening,useInelastic)
            % Initializes the momentum grid, make sure the physicalParamas
            % are properly set and also the TimeGrid tMax
            % It also uses the yMax from MomentumGrid to as a minimum yMax
            % Also make sure relevant grid mode is set in momentumGrid
            % Input:
            % NxiScalingFactor    - scales the predicted Nxi value from simple theory
            %                       considerations; 1 often seems suitable
            % Nxi_min             - minimum possible value for Nxi
            % pMax_default        - the value of pMax if no runaway acceleration occurs
            % pMax_ceiling        - maximum p/mc that we allow 
            % pMaxIncreaseFactor  - uniformly scales the increase in pMax; set to 1 to get
            %                       exactly the predicted value (some margin recommended, ~2)
            % minPMaxMarginFactor - never let pMax be below minPMaxMarginFactor*p0
            %                       (must at least cover the entire distribution, ~1.3)
            if nargin < 2
                useScreening = false;
                useInelastic = false;
            elseif nargin < 3
               useInelastic = false;
            end
            deltaBar = sqrt(min(this.physicalParams.T))/this.reference.TRef;
            tMax = this.timeGrid.tMax/deltaBar^3;
            maxAbsEEc = max(abs(this.physicalParams.EOverEc));
            pc = 1/sqrt(maxAbsEEc-1);
            neTotalOverneFree = this.physicalParams.neTotalOverneFree(1);
            p0 = this.momentumGrid.yMax*this.reference.deltaRef;
            pMax_default = p0;
            if maxAbsEEc > 1
                if 1.5*pc <= pMax_default
                    nuCOvernuCODE = max(this.physicalParams.lnLambdas)/this.reference.lnLambdaRef * 3*sqrt(pi)*this.reference.deltaRef^3/4 *this.physicalParams.neTotalOverneFree(1)* max(this.physicalParams.nBars); % lnLambda variation here too?

                    options = optimoptions('fsolve','Display','off');
                    if useInelastic
                        % Since the logarithmic term in nu_s varies a lot with useInelastic, we need a larger
                        % margin with the predicted pMax, or it breaks down at lowish energy

                        newPMaxPredict = real(fsolve(@(p) p - p0 - (maxAbsEEc/neTotalOverneFree-1)*tMax * nuCOvernuCODE  ...
                                    - pc/2*log( (p/pc +1)./(p/pc-1).*(p0/pc - 1)./(p0/pc+1) ), ...
                                    p0 + (maxAbsEEc/neTotalOverneFree-1)*tMax*nuCOvernuCODE,options));
                        newPMax = max(this.minPMaxMarginFactor*p0, p0 + 1.5*this.pMaxIncreaseFactor*(newPMaxPredict-p0)); % never cut off the solution below p0
                    else
                        newPMaxPredict = real(fsolve(@(p) p - p0 - (maxAbsEEc-1)*tMax * nuCOvernuCODE  ...
                                   - pc/2*log( (p/pc +1)./(p/pc-1).*(p0/pc - 1)./(p0/pc+1) ), ...
                                   p0 + (maxAbsEEc-1)*tMax*nuCOvernuCODE,options)); % 30% above the theoretical momentum a particle with initial momentum p0 will
                                                                                  % reach after being accelerated in the maximum E-field for a duration tMax
                        newPMax = max(this.minPMaxMarginFactor*p0, p0 + this.pMaxIncreaseFactor*(newPMaxPredict-p0)); % never cut off the solution below p0
                    end
                                                                                  % (the formula follows from dp/dt = eEc( E/Ec - 1/v^2), integrated from p0 to p)
                else
                    newPMax = pMax_default;
                end
            else
                newPMax = pMax_default;
            end

            pMax = min(newPMax, this.pMax_ceiling);

            if useScreening 
                ZAtomicNumber = this.physicalParams.species.ZAtomicNumber;
                Ztot = (this.physicalParams.species.nj' * ZAtomicNumber.^2)./ (this.physicalParams.n') ; 
                Zeff = this.physicalParams.Z(:);
                nuDPredict = min(1 + Zeff + 0.3 * (Ztot-Zeff)); % the 0.3 comes from theory, the 5 because the theory doesn't seem to work well? annoying.
                NxiReqPredict = ceil(1.5*this.NxiScalingFactor*sqrt((1+maxAbsEEc ) * pMax/nuDPredict));
            else
                NxiReqPredict =  ceil(1.5*this.NxiScalingFactor* sqrt( (maxAbsEEc+1) * pMax/(min(this.physicalParams.Z)+1) ) );
            end
            Nxi  = max(NxiReqPredict, this.Nxi_min);
            Nxi  = min(Nxi, this.Nxi_max);

            if this.momentumGrid.yGridMode == 6
                yMax      = pMax/this.reference.deltaRef;
                DeltaY0   = deltaBar/(10 + tMax^(1/4)) ;  %in CODE collision times
                DeltaYmax = this.dyTailScalingFactor * DeltaY0*yMax/2 ; % assume that the characteristic momentum scale 
                DeltaY0   = this.dyBulkScalingFactor * DeltaY0;         % of the bulk is y=1 and at yMax is y=yMax/2  

                s0 = this.percentBulk;
                y0 = this.pSwitch;
                if y0 == 0
                    %if gridStepPosition = 0, set to near critical momentum
                    if maxAbsEEc > 1
                        y0 = 0.9*this.reference.deltaRef^(-1)*pc;
                    else
                        y0 = 5;
                    end
                end
                eta = (1-(y0/yMax)^(this.r-1))/(this.r-1);
    %             k   = eta*s0/(1-s0);

                s0_min = fzero(@(s) (this.r*eta^2 + 4*eta+6)*s^3+(this.r*eta^2-4*eta-18)*s^2+18*s-6,0.5);
                if this.r>= 2/3
                    k_crit = 2/this.r*(1+(3*this.r/2-1)^(1/3));
                else
                    k_crit = 2/this.r*(1-sqrt(1-3*this.r/2));
                end
                s0_crit = k_crit/(eta+k_crit);
                if (s0_min > 0) && (s0_min < 0.99*s0_crit) 
                    s0_max = s0_min;
                else
                    s0_max = 0.99*s0_crit;
                end

                if s0 > s0_max
                    % if s0 exceeds s0_min or s0_crit, 
                    % set equal to the lower one
                    s0 = s0_max;
                end


                kappa = (y0/yMax)^this.r;
                rho = DeltaY0/DeltaYmax;

                s0New = s0;
                if rho/kappa > sqrt(6*this.r)-2 %exists a real solution to kOpt below, then overwrite s0 with this choice
                    kOpt = 1/this.r*(2+rho/kappa*(1 - sqrt((2*kappa/rho + 1)^2 - 6*r*kappa^2/rho^2 ) ));
                    sOpt = kOpt/(kOpt+eta);
                    if sOpt < s0_max
                        s0New = sOpt;
                    end

                    % With the below we always spend at least 10% of grid points on the bulk
                    if sOpt < 0.1
                        s0New = 0.1;
                    end
                end

                Ny1 = y0/DeltaY0 * eta/(1-s0New)*(3/eta*(1-s0New)/s0New - 2 + this.r*eta/2*s0New/(1-s0New));
                Ny2 = yMax/DeltaYmax * (yMax/y0)^(this.r-1)*eta/(1-s0New);

                NyNew = round(max(Ny1,Ny2));
            end
            this.momentumGrid.setResolution('yMax',pMax/this.reference.TRef,'Nxi',Nxi,'Ny',NyNew,'gridStepWidth',s0New);
        end
        
    end

end
