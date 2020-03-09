classdef Output < handle
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Settings
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % saveDist - switch if distribution is saved, default true
    properties
        saveDist = true
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Saved Values
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % distrubitions - cell with saved disitributions
    % timesteps - timesteps saved
    % 
    properties (SetAccess = protected)
        distributions = {};
        yc = {};
        nr = {};
        averageEnergyMeV = {};
        averageFastEnergyMeV = {};
        averageREEnergyMeV = {};
        currentDensity = {};
        currentDensityRE = {};
        fracRE = {};
        fracRECurrent = {};
        fracREEnergy = {};
        growthRate = {};
        growthRatePerSecond = {};
        times = {};
        totalEnergyMeV = {};
        totalEnergyJ = {};
        dJdtSI = {};
        dJ_rundtSI = {};
    end

    properties (SetAccess = protected)
        T = {};
        n = {};
        Z = {};
        E = {};
        B = {};
        species = {};
        neTotalOverneFree = {};
        nuees = {};
        deltas = {};
        EOverEc = {};
        EOverED = {};
        EHats = {};
        BHatRef = {};
        nueeBars = {};
        nBars = {};
        veBars = {};
        veBars2 = {};
        veBars3 = {};
        lnLambdas = {};
    end
    methods
        function this = Output()
            
        end

        function save(this, timeIndex, state, f, fbefore)
            %% TODO controll these calulations
            %{
            restMassEnergy = 0.510998928; %MeV
            eG = 4.80320425e-10;
            cG = 2.9979e10; %Gaussian units
            currentConversionFactor = -4/(3*sqrt(pi))*eG*1e-6*state.reference.nRef*state.reference.deltaRef*cG/2.9979e5; %in ampere
            energyMoment = state.momentumGrid.y2.*(state.momentumGrid.gamma-1);
            energyConversionFactorMeV = 4/sqrt(pi)*restMassEnergy*state.reference.nRef;
            totalEnergy = (state.momentumGrid.yWeights.*energyMoment) * f(1:state.momentumGrid.Ny);
            this.yc{end+1} = 1/sqrt(abs(state.physicalParams.EOverEc(timeIndex))-1);
            this.nr{end+1} = (state.momentumGrid.yWeights.*state.momentumGrid.y2.*(state.momentumGrid.y>this.yc{end}))*f(1:state.momentumGrid.Ny)*state.reference.nRef;
            
            this.averageEnergyMeV{end+1} = energyConversionFactorMeV*totalEnergy/state.physicalParams.n(timeIndex);
            %averageFastEnergyMeV = (state.momentumGrid.yWeights.*energyMoment.*(state.y>this.yc{end})) * f(1:state.momentumGrid.Ny); % this is runaway current right now, not fast particles
            this.averageREEnergyMeV{end+1} = energyConversionFactorMeV*(state.momentumGrid.yWeights.*energyMoment.*(state.momentumGrid.y>this.yc{end})) * f(1:state.momentumGrid.Ny)/this.nr{end};

            this.currentDensity{end+1} = currentConversionFactor*(state.momentumGrid.xy2.*state.momentumGrid.yWeights)*f(state.momentumGrid.Ny+(1:state.momentumGrid.Ny));
            this.currentDensityRE{end+1} = currentConversionFactor*(state.momentumGrid.xy2.*state.momentumGrid.yWeights.*(state.momentumGrid.y>this.yc{end}))*f(state.momentumGrid.Ny+(1:state.momentumGrid.Ny));
            this.fracRE{end+1} = this.nr{end}/state.physicalParams.n(timeIndex);
            this.fracRECurrent{end+1} = this.currentDensityRE{end}/this.currentDensity{end};
            this.fracREEnergy{end+1} = this.averageREEnergyMeV{end}*this.nr{end}/this.averageEnergyMeV{end}/state.physicalParams.n(timeIndex);

            this.times{end+1} = state.timeGrid.timesteps(timeIndex);
            this.totalEnergyMeV{end+1} = this.averageEnergyMeV{end}*state.physicalParams.n(timeIndex);
            this.totalEnergyJ{end+1} = this.totalEnergyMeV{end}*1.6021773e-13;
            if timeIndex > 1 
                this.growthRate{end+1} = (this.nr{end}-...
                    (state.momentumGrid.yWeights.*state.momentumGrid.y2.*...
                    (state.momentumGrid.y>1/sqrt(abs(state.physicalParams.EOverEc(timeIndex-1))-1)))*...
                    fbefore(1:state.momentumGrid.Ny)...
                    *state.reference.nRef)/state.timeGrid.dts(timeIndex-1);
                this.growthRatePerSecond{end+1}= this.growthRate{end}*state.reference.nueeRef;
                this.dJdtSI{end+1} = (this.currentDensity{end}-...
                    currentConversionFactor*(state.momentumGrid.xy2.*state.momentumGrid.yWeights)*...
                    fbefore(state.momentumGrid.Ny+(1:state.momentumGrid.Ny)))/state.timeGrid.dts(timeIndex-1)*state.reference.nueeRef;
                this.dJ_rundtSI{end+1}= (this.currentDensityRE{end} - ...
                    currentConversionFactor*(state.momentumGrid.xy2.*state.momentumGrid.yWeights.*...
                    (state.momentumGrid.y>1/sqrt(abs(state.physicalParams.EOverEc(timeIndex-1))-1)))*...
                    fbefore(state.momentumGrid.Ny+(1:state.momentumGrid.Ny)))/state.timeGrid.dts(timeIndex-1)*state.reference.nueeRef;
            else
                this.growthRate{end+1}= 0;
                this.growthRatePerSecond{end+1} = 0;
                this.dJdtSI{end+1} = 0;
                this.dJ_rundtSI{end+1} = 0;
            end
            %}
            %%
            this.T{end+1} = state.physicalParams.T(timeIndex);
            this.n{end+1} = state.physicalParams.n(timeIndex);
            this.Z{end+1} = state.physicalParams.Z(timeIndex);
            this.E{end+1} = state.physicalParams.E(timeIndex);
            this.B{end+1} = state.physicalParams.B;
            this.species{end+1} = state.physicalParams.species;
            this.neTotalOverneFree{end+1} = state.physicalParams.neTotalOverneFree(timeIndex);
            this.nuees{end+1} = state.physicalParams.nuees(timeIndex);
            this.deltas{end+1} = state.physicalParams.deltas(timeIndex);
            this.EOverEc{end+1} = state.physicalParams.EOverEc(timeIndex);
            this.EOverED{end+1} = state.physicalParams.EOverED(timeIndex);
            this.EHats{end+1} = state.physicalParams.EHats(timeIndex);
            this.BHatRef{end+1} = state.physicalParams.BHatRef;
            this.nueeBars{end+1} = state.physicalParams.nueeBars(timeIndex);
            this.nBars{end+1} = state.physicalParams.nBars(timeIndex);
            this.veBars{end+1} = state.physicalParams.veBars(timeIndex);
            this.veBars{end+1} = state.physicalParams.veBars(timeIndex);
            this.veBars{end+1} = state.physicalParams.veBars(timeIndex);
            this.lnLambdas{end+1} = state.physicalParams.lnLambdas(timeIndex);
            
            if this.saveDist
                this.addDistribution(f,state,timeIndex);
            end
        end


        function addDistribution(this,f,state,saveIndex)
            % STATE adds a distribution `f` at index `saveIndex`.
            this.distributions{end+1} = Distribution(f, ...
                            state,...
                            saveIndex);
        end



        function [f, times] = getDistributions(this)
         % Returns all distributions `f` at all times `times`. The
         % output is sorted after times, and normalization is broken.
         % Distribution is given as (Ns)^-3 m^-3 and integrates to the density n 
         % when integreted over all momentum
            N = length(this.distributions);
            times = zeros(1,N);
            for i = 1:N
                times(i) = this.distributions{i}.time;
            end
            [times,ts] = sort(times);
            I = 1:N;
            I = I(ts);
            f = cell(1,N);
            for i = 1:N
                m0 = 9.10938356e-31; %SI
                c = 299792458; %SI
                f{i} = this.distributions{I(i)}.f*this.distributions{I(i)}.nRef...
                    /(sqrt(pi)*m0*this.distributions{I(i)}.deltaRef*c)^3;
            end
        end


        function [p, times] = getMomentumVectors(this)
            % Returns all momentum vectors `p` at all times `times`. The
            % output is sorted after times, and normalization is broken.
            % Momentum is given as Ns.
            N = length(this.distributions);
            times = zeros(1,N);
            for i = 1:N
                times(i) = this.distributions{i}.time;
            end
            [times,ts] = sort(times);
            I = 1:N;
            I = I(ts);
            p = cell(1,N);
            for i = 1:N
                m0 = 9.10938356e-31; %SI
                c = 299792458; %SI
                p{i} = this.distributions{I(i)}.momentumGrid.y*this.distributions{I(i)}.deltaRef*c*m0;
            end
        end
    end

end
