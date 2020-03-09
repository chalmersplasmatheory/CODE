classdef Distribution < handle
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    % Distribution
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % f 
    % momentumGrid
    % time
    
    %Design note, momentum grid is not checked for changes when saved to
    %this class. Therefore if things change in memory to which this points
    %to, it will be unnoticed and wrong.?
    properties (SetAccess = immutable)
        f
        momentumGrid
        time
        % for convinient plot scripts to only take distibution
        T
        n
        Z
        E
        B
        species
        neTotalOverneFree
        nuees
        deltas
        EOverEc
        EOverED 
        EHats 
        BHatRef 
        nueeBars 
        nBars 
        veBars 
        veBars2 
        veBars3 
        lnLambdas 
        
        %These reference values should be replaced with reference, or just
        %save the distribution as unnormalized data?
        TRef
        nRef
        deltaRef
        lnLambdaRef
        nueeRef
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%  Constructor      %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
       
        function this = Distribution(f,state,timeIndex)
            this.f = f;
            this.momentumGrid = state.momentumGrid;
            this.time = state.timeGrid.timesteps(timeIndex);
            
            %Reference
            this.TRef = state.momentumGrid.reference.TRef;
            this.nRef = state.momentumGrid.reference.nRef;
            [this.deltaRef,this.lnLambdaRef,this.nueeRef,~] = getDerivedParameters(this.TRef,this.nRef,0);
            
            
            this.T = state.physicalParams.T(timeIndex);
            this.n = state.physicalParams.n(timeIndex);
            this.Z = state.physicalParams.Z(timeIndex);
            this.E = state.physicalParams.E(timeIndex);
            this.B = state.physicalParams.B;
            this.species = state.physicalParams.species;
            this.neTotalOverneFree = state.physicalParams.neTotalOverneFree(timeIndex);
            this.nuees = state.physicalParams.nuees(timeIndex);
            this.deltas = state.physicalParams.deltas(timeIndex);
            this.EOverEc = state.physicalParams.EOverEc(timeIndex);
            this.EOverED = state.physicalParams.EOverED(timeIndex);
            this.EHats = state.physicalParams.EHats(timeIndex);
            this.BHatRef = state.physicalParams.BHatRef;
            this.nueeBars = state.physicalParams.nueeBars(timeIndex);
            this.nBars = state.physicalParams.nBars(timeIndex);
            this.veBars = state.physicalParams.veBars(timeIndex);
            this.veBars = state.physicalParams.veBars(timeIndex);
            this.veBars = state.physicalParams.veBars(timeIndex);
            this.lnLambdas = state.physicalParams.lnLambdas(timeIndex);
        end
        
    end
    
end
