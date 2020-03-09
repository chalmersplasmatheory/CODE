classdef TimeGrid < handle
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Time advance:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dt - timestep size (normalized)
    % tMax - maximum time (minimum is so far always 0) (normalized)
    % timeStepMode - specifies how the time step vector is generated
    %       0 = Constant time step (dt), in units of timeUnit
    %       2 = Logarithmic -- use progressively longer time steps. Usefull
    %           for convergence towards a steady state. dt is used for the
    %           first time step.
    %       3 = Stepwise logarithmic -- like 2, except that several time
    %           steps are taken with each step length. This is to avoid
    %           rebuilding the matrix in every time step.
    %       4 = fully user supplied time vector
    % logGridScaling - step size scaling for the logarithmic grid
    % logGridSubSteps - how many steps to take for each time step length
    %                   with timeStepMode 3
    % logGridMaxStep - maximum time step allowed in timeStepModes 2 and 3
    % timesteps - actual times the distribution is calculated in (normalized)
    % dts - vector of all small timechanges in (normalized)
    % dtsHasChanged - vector contantaining if the a element in dts vector
    %                   has changed or not
    % nTimeSteps - number of times where the parameters are defined
    %               (actually one greater than the number of timesteps
    %               that should be taken since we get the initial timestep for free)
    properties (SetAccess = protected)
        %Time advance
        dt = 1; %Distance between different time steps in uniform grid. In case of nonuniform, determines the order of magnitude of timedifference
        tMax = 1;% Maximum time to simulate to. This will always be exact (to machine precision)
        timeStepMode = 0; %Switch for time spacing
        logGridScaling = sqrt(2); % Scaling in logarithmic time steps
        logGridSubSteps = 5; %In stepwise logarithmic grid mode, how many uniform step per increase
        logGridMaxStep = 100000; % how many steps are tolerated
        timesteps %Vector containing all the different times
        dts %Vector containing all the dts used, namely diff(timesteps)
        dtsHasChanged %vector containing if dts has changed
        nTimeSteps %number of timeSteps
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%
    % Misc
    % %%%%%%%%%%%%%%%%%%%%%%%%
    %Design note: intention is that physicalParameters containing an instance of
    % this class also is contained in this object, think of it as pairs. 
    % Right now it is possible to 'hack' the construction by creating a
    % TimeGrid object, and a PhysicalParams. Then setting the TimeGrid in
    % the PhysicalParams. Afterwards it is possible to create a new 
    % PhysicalParams and set the already
    % created TimeGrid in the new PhysicalParams. The firstly created PhysicalParams 
    % now has a TimeGrid object in it with a PhysicalParams object in it which 
    % is not pointing to itself. The 'pair' structure is then broken. Two
    % fixes are availible: one seperating the TimeGrid to a copy of itself
    % (new pointer) and using the copy for the old TimeGrid or save both
    % TimeGrid objects in the same reference object.
    properties (SetAccess = protected)
        physicalParams 
        reference
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%    Constructor    %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function this = TimeGrid(reference,varargin)
            for k = 1:2:(nargin-1)
                if isa(this.(varargin{k}),'Reference')
                    error('This function does not yet support an update of the reference object in other than as first argument.')
                elseif isa(this.(varargin{k}),'PhysicalParams')
                    this.setPhysicalParams(this,varargin{k+1})
                else
                    this.(varargin{k}) = varargin{k+1};
                end
            end
            this.reference = reference;
            reference.setTimeGrid(this);
            
            this.initializeTimeGrid();
            if ~isempty(this.physicalParams)
                this.physicalParams.interpolatePhysicalParams();
            end
        end

    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%  physicalParams   %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function setPhysicalParams(this,phP)
            if (isa(phP,'PhysicalParams'))
                if(phP.timeGrid == this) %check if timeGrid in phP points to this object
                    this.physicalParams = phP;
                else
                    error('PhysicalParams input should contain a TimeGrid pointer to this instance of TimeGrid.')
                end
            else
                error('Input must be PhysicalParams.')
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%    time grid     %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%       (set)      %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function setResolution(this,varargin)
            for k = 1:2:(nargin-1)
                if isa(this.(varargin{k}),'Reference')
                    error('This function does not yet support an update of the reference object.')
                elseif isa(this.(varargin{k}),'PhysicalParams')
                    this.setPhysicalParams(this,varargin{k+1})
                else
                    this.(varargin{k}) = varargin{k+1};
                end
            end
            this.initializeTimeGrid();
            if ~isempty(this.physicalParams)
                this.physicalParams.interpolatePhysicalParams();
            end
        end
        
        function setdt(this, dt)
            this.dt = dt;
            this.initializeTimeGrid();
            if ~isempty(this.physicalParams)
                this.physicalParams.interpolatePhysicalParams();
            end
        end

        function settMax(this, tMax)
            this.tMax = tMax;
            this.initializeTimeGrid();
            if ~isempty(this.physicalParams)
                this.physicalParams.interpolatePhysicalParams();
            end
        end

        function settimeStepMode(this, timeStepMode)
            this.timeStepMode  = timeStepMode;
            this.initializeTimeGrid();
            if ~isempty(this.physicalParams)
                this.physicalParams.interpolatePhysicalParams();
            end
        end

        function setlogGridScaling(this, logGridScaling)
            this.logGridScaling  = logGridScaling;
            this.initializeTimeGrid();
            if ~isempty(this.physicalParams)
                this.physicalParams.interpolatePhysicalParams();
            end
        end

        function setlogGridSubSteps(this, logGridSubSteps)
            this.logGridSubSteps = logGridSubSteps;
            this.initializeTimeGrid();
            this.physicalParams.interpolatePhysicalParams();
        end

        function setlogGridMaxStep(this, logGridMaxStep)
            this.logGridMaxStep  = logGridMaxStep;
            this.initializeTimeGrid();
            if ~isempty(this.physicalParams)
                this.physicalParams.interpolatePhysicalParams();
            end
        end

        function settimesteps(this, timesteps, timeUnit)
            switch timeUnit
                case 's' 
                    timesteps = timesteps*this.reference.nueeRef;
                case 'ms'
                    timesteps = timesteps*this.reference.nueeRef*1e-3;
                case 'normalized'
                otherwise
                    error('Unknown time unit supplied')
            end
            this.dt = NaN; %to make sure that we get errors if dt is used, should use dts instead
            this.tMax = this.timesteps(end);
            this.timeStepMode = 4; %Switch for time spacing
            this.timesteps = timesteps; %Vector containing all the different times
            this.dts = diff(timesteps); %Vector containing all the dts used, namely diff(timesteps)
            this.dtsHasChanged = [1 abs(diff(this.dts))>1e-13]; %vector containing if dts has changed
            this.nTimeSteps = numel(timesteps)-1; %number of timeSteps
        end

        % here f can be used so that it is saved at only specific times, maybe in initializeGrid (todo)
        function initializeTimeGrid(this)
            if this.timeStepMode == 4
               warning('Note that if you supply your own timesteps, other changes to timeGrid will not change')
               return
            end
            this.timesteps = [];
            this.dts = [];
            switch this.timeStepMode
                case 0 %Standard, linear in reference collision times
                    this.timesteps = 0:this.dt:this.tMax;
                    this.dts = this.dt*ones(size(this.timesteps));
                case 2 %Logarithmic
                    incFact = this.logGridScaling;%default sqrt(2); %Determine n steps
                    nSteps = log((incFact-1)*this.tMax/this.dt+1)/log(incFact);
                    this.timesteps = logspace(log10(this.dt),log10(this.tMax),ceil(nSteps)+1);

                    %Keep the time steps length below some value
                    id1 = find(this.timesteps>this.logGridMaxStep,1);
                    if ~isempty(id1)
                        this.timesteps = [this.timesteps(1:id1-1),linspace(this.timesteps(id1),this.timesteps(end),ceil((this.timesteps(end)-this.timesteps(id1))/this.logGridMaxStep))];
                    end
                    this.timesteps = [0 this.timesteps];
                    this.dts = diff(this.timesteps);
                case 3 %Logarithmic, step-wise
                    incFact = this.logGridScaling;%sqrt(2); %Determine n logarithmic steps
                    DT = this.dt;
                    nSubSteps = this.logGridSubSteps;%5;
                    this.timesteps = zeros(1,ceil(this.tMax/DT)+1);%the maximum possible length
                    iSteps = 0;
                    remainingSteps = this.tMax/DT;
                    while remainingSteps>nSubSteps
                        ind = (1:nSubSteps) + (iSteps)*nSubSteps+1;
                        this.timesteps(ind) = DT:DT:DT*nSubSteps;
                        DT = DT*incFact;
                        this.timesteps(ind) = this.timesteps(ind(1)-1)+(DT:DT:DT*nSubSteps);
                        remainingSteps = (this.tMax-this.timesteps(ind(end)))/(DT*incFact);
                        iSteps = iSteps+1;
                    end
                    %last: make sure we end at tMax.
                    nRemainingSteps = ceil(remainingSteps);
                    if ~isempty(ind)
                        DT = (this.tMax-this.timesteps(ind(end)))/nRemainingSteps;
                    else
                        DT = this.tMax/nRemainingSteps;
                    end
                    ind = (1:nRemainingSteps) + (iSteps)*nSubSteps+1;
                    this.timesteps(ind) = this.timesteps(ind(1)-1)+(DT:DT:DT*nRemainingSteps);
                    this.timesteps = this.timesteps(1:ind(end));
                    %Keep the time steps length below some value
                    id1 = find(this.timesteps>this.logGridMaxStep,1);
                    if ~isempty(id1)
                        this.timesteps = [this.timesteps(1:id1-1),linspace(this.timesteps(id1),this.timesteps(end),ceil((this.timesteps(end)-this.timesteps(id1))/this.logGridMaxStep))];
                    end
                    
                    this.dts = diff(this.timesteps);
                case 4 %self supplied. nothing todo
                otherwise
                    error('Invalid time step mode');
            end
            this.dtsHasChanged = abs(diff(this.dts)) > 1e-13;
            this.dtsHasChanged = [1,this.dtsHasChanged]; %diff removes one element, first is same as initilizing process
            this.nTimeSteps = numel(this.timesteps);         
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%    time grid     %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%       (get)      %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		function dt = getdt(this,timeUnit)
            dt = this.dt;
            if nargin >1
               switch timeUnit
                   case 's'
                       dt = dt/this.nueeRef;
                   case 'ms'
                       dt = dt*1e-3/this.nueeRef;
                   case 'normalized'
                   otherwise
                       error('Unknown time unit.')
               end
            end
        end

        function tMax = gettMax(this,timeUnit)
			tMax = this.tMax;
            if nargin >1
               switch timeUnit
                   case 's'
                       tMax = tMax/this.nueeRef;
                   case 'ms'
                       tMax = tMax*1e-3/this.nueeRef;
                   case 'normalized'
                   otherwise
                       error('Unknown time unit.')
               end
            end
        end

		function timeStepMode  = gettimeStepMode(this)
			timeStepMode  = this.timeStepMode;
		end

		function logGridScaling  = getlogGridScaling(this)
			logGridScaling  = this.logGridScaling;
		end

		function logGridSubSteps  = getlogGridSubSteps(this)
			logGridSubSteps  = this.logGridSubSteps;
		end

		function logGridMaxStep  = getlogGridMaxStep(this)
			logGridMaxStep  = this.logGridMaxStep;
        end

		function timesteps = gettimesteps(this,timeUnit)
			timesteps = this.timesteps;
            if nargin >1
               switch timeUnit
                   case 's'
                       timesteps = timesteps/this.nueeRef;
                   case 'ms'
                       timesteps = timesteps*1e-3/this.nueeRef;
                   case 'normalized'
                   otherwise
                       error('Unknown time unit.')
               end
            end
        end

        function dts = getdts(this,pts,timeUnit)
            if nargin == 1
                dts = this.dts;
            else
                dts = this.dts(pts);
            end
            if nargin >2
               switch timeUnit
                   case 's'
                       dts = dts/this.nueeRef;
                   case 'ms'
                       dts = dts*1e-3/this.nueeRef;
                   case 'normalized'
                   otherwise
                       error('Unknown time unit.')
               end
            end
        end

        function dtsHasChanged = getdtsHasChanged(this)
			dtsHasChanged = this.dtsHasChanged;
        end

        function nTimeSteps = getnTimeSteps(this)
			nTimeSteps = this.nTimeSteps;
		end
    end
    
end
