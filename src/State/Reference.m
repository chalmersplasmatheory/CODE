classdef Reference < handle
    %REFERENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    %  Reference values
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    % TRef - Reference temperature
    % nRef - Reference density
    % nueeRef - Reference collision frequency
    % deltaRef - Reference thermal velocity over speed of light
    % lnLambdaRef - Reference coulumb logarithm    
    
    properties (SetAccess = protected)
        TRef, nRef, nueeRef, deltaRef, lnLambdaRef 
    end
    
    %Design note: intention is that all objects containing an instance of
    % this class also is contatined in this object, think of it as pairs. 
    % Question is if more than
    % one instance of each class should be able to share a reference object.
    % Right now it is possible to 'hack' the construction by creating a
    % Reference object, and for example a TimeGrid. Then setting the reference in
    % the TimeGrid. Afterwards it is possible to create a new TimeGrid and set the already
    % created Reference in the new TimeGrid. The firstly created TimeGrid 
    % now has a Reference object in it with a TimeGrid object in it which 
    % is not pointing to itself. The 'pair' structure is then broken. Two
    % fixes are availible: one seperating the reference to a copy of itself
    % (new pointer) and using the copy for the old TimeGrid or save both
    % TimeGrid objects in the same reference object.
    properties (SetAccess = protected)
        physicalParams
        momentumGrid
        timeGrid
    end
    
    methods
        
        function this = Reference(TRef,nRef)
            this.TRef = TRef;
            this.nRef = nRef;
            [this.deltaRef,this.lnLambdaRef,this.nueeRef,~] = getDerivedParameters(this.TRef,this.nRef,0);
        end
        
        function setPhysicalParams(this,phP)
            if (isa(phP,'PhysicalParams'))
                if(phP.reference == this) %check if Reference in phP points to this object
                    this.physicalParams = phP;
                else
                    error('PhysicalParams input should contain a Reference object pointer to this instance of Reference.')
                end
            else
                error('Input must be PhysicalParams.')
            end
        end
        
        function setMomentumGrid(this,mg)
            if (isa(mg,'MomentumGrid'))
                if (mg.reference == this)
                    this.momentumGrid = mg;
                else
                    error('MomentumGrid must contain this instance of the Reference class.')
                end
            else
                error('Input must be a MomentumGrid.')
            end
        end
        
        function setTimeGrid(this,tg)
            if isa(tg, 'TimeGrid')
                if (tg.reference == this)
                    this.timeGrid = tg;
                else
                    error('TimeGrid must contain this instance of the Reference class.')
                end
            else
                error('Input must be a TimeGrid.')
            end
        end
        
    end
    
    methods
       
        function updateReferenceVals(this,TRef,nRef)
            %Updates reference values and renormalizes values in owned
            %TimeGrid, MomentumGrid and PhysicalParams
            [deltaRefNew,lnLambdaRefNew,nueeRefNew,~] = getDerivedParameters(TRef,nRef,0);
            convFac = 1/this.nueeRef*nueeRefNew; %Time convertion factor
            %Update Time Grid
            if ~isempty(this.timeGrid)
                if this.timeGrid.timeStepMode == 4 %only use set timesteps if we have used it before
                    this.timeGrid.settimesteps(this.timeGrid.timesteps*convFac,'normalized');
                else
                    this.timeGrid.setdt(this.timeGrid.dt*convFac)
                    this.timeGrid.settMax(this.timeGrid.tMax*convFac)
                end
            end
            %Update PhysicalParams raw data !!BEFORE ReferenceVals Are Set!
            if ~isempty(this.physicalParams)
                this.physicalParams.updateRawTimes(nueeRefNew);
            end
            
            this.TRef = TRef;
            this.nRef = nRef;
            this.deltaRef = deltaRefNew;
            this.lnLambdaRef = lnLambdaRefNew;
            this.nueeRef = nueeRefNew;
            
            %Update momentum grid !!AFTER new values
            if ~isempty(this.momentumGrid)
                this.momentumGrid.initializeMomentumGrid();
            end
            %Update PhysicalParams dependent params !!After new reference
            %are set
            if ~isempty(this.physicalParams)
                this.physicalParams.calcDepParams();
            end
        end
        
    end
end

