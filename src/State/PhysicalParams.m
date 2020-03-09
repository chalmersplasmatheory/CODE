classdef PhysicalParams <  handle
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%                   Description of properties                     %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % %%%%%%%%%%%%%%%%%%%%%%
    %  Physical quantities:
    % %%%%%%%%%%%%%%%%%%%%%%
    %
    % T - Plasma temperature (eV), given in timesteps
    % n - Plasma density (m^{-3}), given in timesteps
    % Z - Plasma effective charge, given in timesteps
    % E - Electric field (V/m), given in timesteps
    % B - (optional) Magnetic field (T), only used for calculating
    %     synchrotron radiation reaction. If a magnetic field strength is
    %     provided, a term describing momentum loss due to synchrotron
    %     emission is included in the kinetic equation. If B is empty it is
    %     not included.
    % timeGrid - object of TimeGrid
    % species - a struct used with the species properties so that
    %     species has the fields
    %     nj, Z0NetCharge, ZAtomicNumber, times:
    %     nj(jspecies, times) - matrix containing number density of each
    %       species as a function of time (m^{-3})
    %     ZAtomicNumber - atomic number of each species. row vector, not a
    %     function of time
    %     Z0NetCharge - net charge for each species. row vector, not a
    %     function of time
    %     times - time steps. If constant, just set i to 1.
    %     NB species cannot be set from SetPhysicalParameters.
    %     NB overrides specified Z and n since nj and Z0NetCharge
    %     specifies the effective charge. (Z is the fully screened value)
    %     NB only used if useScreening ~= 0
    % rawT, rawn, rawZ, rawE, rawB - same as without raw except not changed
    % and given in raw** and is the raw input that the user supplied.
    % rawtT - times where rawT is defined, user supplied (normalized)
    % rawtn - times where rawn is defined, user supplied (normalized)
    % rawtZ - times where rawZ is defined, user supplied (normalized)
    % rawtE - times where rawE is defined, user supplied (normalized)
    % neTotalOverneFree - vector containing fraction density of total
    %       electrons compared to free electrons at all timesteps.
    %       Used when specicies is used (for impurities and boltzmann
    %       operator, and bound electrons are present). If no species are
    %       used or no source is used, this is set to 1.
    % nuees - Collisional frequency of electron electron collisions
    % deltas -  thermal speed over speed of light
    % EOverEc - Electric field over the Critical Electrical field
    % EOverED - Electric field over the Dreicer Electrical field
    % lnLambdas - Coloumb logarithm
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Miscellaneous:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % interpolationMethod - Interpolation method to interpolate between
    %                       given values and values at timestepping times.
    %                       Used in interp1 (default previous) NOTE: SINCE
    %                       NEAREST IS DEFAULT, T WILL DROP AT HALF THE
    %                       TIME YOU SPECIFY WHERE T SHOULD DROP


    properties (SetAccess = protected)
        %Physical
        T, n, Z = 0, E = 0, B = 0
        species=[];
        
        % We also save the raw input to use if reinterpolation is done,
        % these are only used internally
        rawT = 0, rawn = 0, rawZ = 0,rawE = 0
        rawtT = 0, rawtn = 0, rawtZ = 0, rawtE = 0 % normalized
        rawspecies
        
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
        
        %misc
        interpolationMethod = 'previous';
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%
    % Time Grid
    % %%%%%%%%%%%%%%%%%%%%%%%%
    %immutable so that timeGrid can be sure that if passed PhysicalParams
    %to TimeGrid contains reference to it self, it is a valid
    %physicalParams input
    properties (SetAccess = protected) 
        timeGrid 
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%
    %  Reference values
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    % TRef - Reference temperature
    % nRef - Reference density
    % nueeRef - Reference collision frequency
    % deltaRef - Reference velocity over speed of light
    % lnLambdaRef - Reference coulumb logarithm    
    % 
    % all contained in a single object `reference`
    properties (SetAccess = protected)
        reference
    end

    properties (Constant)
        VERSION = 1.0;
    end

    methods
        
        function this = PhysicalParams(reference,timeGrid,varargin)
        % PHYSICALPARAMS initiates a PhysicalParams by setting values to
        % to what is specified. First two arguments are Reference object, 
        % and TimeGrid. The rest are optional and follow standard matlab
        % syntax (value, 'name')
            if isa(timeGrid,'TimeGrid') 
                this.timeGrid = timeGrid; 
            end
            if isa(reference, 'Reference')
                this.reference = reference;
            end
            for k = 1:2:(nargin-2)
                if isa(this.(varargin{k}),'Reference')
                    error('This function does not yet support an update of the reference object in other than as first argument.')
                elseif isa(this.(varargin{k}),'TimeGrid')
                    error('This function does not yet support an update of the TimeGrid object in other than as second argument.')
                else
                    this.(varargin{k}) = varargin{k+1};
                end
            end
            timeGrid.setPhysicalParams(this);
            reference.setPhysicalParams(this);
            this.interpolatePhysicalParams;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                        %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Physical properties    %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%       (set)            %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %(todo) user thinks it should supply 'T','n','Z' ... as variable
        %names but should in contrast supply 'rawT', 'rawn' ... as
        %variables
        function setParams(this,varargin)
            for k = 1:2:(nargin-1)
                if isa(this.(varargin{k}),'Reference')
                    error('This function does not yet support an update of the Reference.')
                elseif isa(this.(varargin{k}),'TimeGrid')
                    error('This function does not yet support an update of the TimeGrid object.')
                else
                    this.(['raw' varargin{k}]) = varargin{k+1};
                end
            end
            this.setB(this.B);
            this.interpolatePhysicalParams();
        end
        
        function setPhysicalParameters(o,T,n,Z,E,varargin)
            % Usage:
            % o.SetPhysicalParameters(T,n,Z,E)
            % o.SetPhysicalParameters(T,n,Z,E,B)
            % o.SetPhysicalParameters(T,n,Z,E,B,tT,tn,tZ,tE)
            o.rawT = T;
            o.rawn = n;
            o.rawZ = Z;
            o.rawE = E;
            o.rawtT = 0;
            o.rawtn = 0;
            o.rawtZ = 0;
            o.rawtE = 0;
            switch nargin
                case 6
                    o.setB(varargin{1});
                case 9
                    o.rawtT = varargin{1};
                    o.rawtn = varargin{2};
                    o.rawtZ = varargin{3};
                    o.rawtE = varargin{4};
                case 10
                    o.setB(varargin{1});
                    o.rawtT = varargin{2};
                    o.rawtn = varargin{3};
                    o.rawtZ = varargin{4};
                    o.rawtE = varargin{5};
            end
            o.interpolatePhysicalParams();
        end
                
        function setspecies(this,species)
            this.rawspecies = species;
            disp('Using species, thus overwriting n, and Z.')
            nj = species.nj;
            Z0NetCharge = species.Z0NetCharge;
            
            %new electron density: sum over all species => a new
            %time-dependent vector
            this.rawn = Z0NetCharge' * nj;
            
            %calculate effective charge
            Zeff = (Z0NetCharge'.^2 * nj)./this.rawn;
            this.rawZ = Zeff;
            
            if length(species.times)>=1
                this.rawtn = species.times;
                this.rawtZ = species.times;
            else
                this.rawtn = 0;
                this.rawtZ = 0;
            end
            this.interpolatePhysicalParams();
        end
        
        function setT(this,T,tT)
            this.rawT = T;
            this.rawtT = 0;
            if nargin > 2
                this.rawtT = tT;
            end
            this.interpolatePhysicalParams();
        end
        
        function setn(this,n,tn)
            if ~isempty(this.species)
                warning('n will not be set since species controlls n. Use setspecies([]) to clear species if no species will be used')
                return
            end
            this.rawn = n;
            this.rawtn = 0;
            if nargin > 2
                this.rawtn = tn;
            end
            this.interpolatePhysicalParams();
        end
        
        function setZ(this,Z,tZ)
            if ~isempty(this.species)
                warning('Z will not be set since species controlls Z. Use setspecies([]) to clear species if no species will be used')
                return
            end
            this.rawZ = Z;
            this.rawtZ = 0;
            if nargin > 2
                this.rawtZ = tZ;
            end
            this.interpolatePhysicalParams();
        end
        
        function setE(this,E,tE)
            this.rawE = E;
            this.rawtE = 0;
            if nargin > 2
                this.rawtE = tE;
            end
            this.interpolatePhysicalParams();
        end
        
        function setB(this,B)
            this.B = B;
            [~,~,~,this.BHatRef] = getDerivedParameters(this.reference.TRef,this.reference.nRef,this.B);
        end
        
        function setinterpolationMethod(this,interpolationMethod)
            this.interpolationMethod = interpolationMethod;
            this.interpolatePhysicalParams();
        end
        
        function calcDepParams(this)
            [this.deltas,this.lnLambdas,this.nuees,~] = getDerivedParameters(this.T,this.n,[]);
            [this.EOverEc,this.EOverED,~] = getNormalizedEFields(this.E,this.T,this.n);
            [EOverEcTRefnRef, ~, ~] = getNormalizedEFields(this.E,this.reference.TRef,this.reference.nRef);
            this.EHats = EOverEcTRefnRef .* this.reference.deltaRef ^ 2 * 3 * sqrt(pi) / 4;
            this.nueeBars = this.nuees/this.reference.nueeRef;
            this.nBars = this.n/this.reference.nRef;
            this.veBars = this.deltas/this.reference.deltaRef;
            this.veBars2 = this.veBars.^2;
            this.veBars3 = this.veBars.^3;
            % if there are partially ionized atoms, include all bound
            % electrons
            if ~isempty(this.species) %% || useScreening %<- A: was this before
                it = 1;
                for nj_t = this.species.nj
                    this.neTotalOverneFree(it) = (this.species.ZAtomicNumber' * nj_t)...
                    /( this.species.Z0NetCharge'*nj_t);
                    it = it+1;
                end
            else
                this.neTotalOverneFree = ones(size(this.T)); %A: we need it to be the size of the number of timesteps, which T,E etc is
            end
        end
        
        function interpolatePhysicalParams(this)
            %Temperature
            times = this.timeGrid.timesteps;
            if length(times)>1
                if times(end) > this.rawtT(end)
                    this.T = interp1([this.rawtT, times(end)],[this.rawT, this.rawT(end)],times,this.interpolationMethod);
                else
                    this.T = interp1(this.rawtT, this.rawT,times,this.interpolationMethod);
                end
                %Electric field
                if times(end) > this.rawtE(end)
                    this.E = interp1([this.rawtE, times(end)],[this.rawE this.rawE(end)],times,this.interpolationMethod);
                else
                    this.E = interp1(this.rawtE, this.rawE,times,this.interpolationMethod);
                end
                %Density
                if times(end) > this.rawtn(end)
                    this.n = interp1([this.rawtn, times(end)],[this.rawn this.rawn(end)],times,this.interpolationMethod);
                else
                    this.n = interp1(this.rawtn, this.rawn,times,this.interpolationMethod);
                end
                %charge
                if times(end) > this.rawtZ(end)
                    this.Z = interp1([this.rawtZ, times(end)],[this.rawZ this.rawZ(end)],times,this.interpolationMethod);
                else
                    this.Z = interp1(this.rawtZ, this.rawZ,times,this.interpolationMethod);
                end
                %Species, done differently, now times need to span the
                %interpolation
                if ~isempty(this.rawspecies)
                    this.species = struct();
                    %nj = zeros(nSpecies,numel(this.grid.timesteps)); previous
                    %used this in loop and then this.plasma.species.nj = nj might be bugging now
                    for iSpecies = 1:size(this.rawspecies.nj,1) %number of species
                        if(any(abs(diff(this.rawspecies.nj(iSpecies,:))) > 1e-13))
                            this.species.nj(iSpecies,:) = interp1(this.rawspecies.times,this.rawspecies.nj(iSpecies,:),times,this.interpolationMethod);
                        else
                            this.species.nj(iSpecies,:) = interp1([times(1) times(end)],this.rawspecies.nj(iSpecies,1)*ones(1,2),times,'nearest');
                        end
                    end
                    this.species.ZAtomicNumber = this.rawspecies.ZAtomicNumber;
                    this.species.Z0NetCharge = this.rawspecies.Z0NetCharge;
                end
            end
            this.calcDepParams();
            if ~isempty(this.species)
                this.species.times = this.timeGrid.timesteps;
            end
        end
        
    end
    
    methods
        %
        function updateRawTimes(this, nueeRefNew)
            %PHYSICALPARAMS DO NOT USE THIS FUNCTION, IT IS ONLY FOR
            %REFERENCE CLASS. NO SOLUTION FOUND WHERE ONLY REFERENCE CLASS
            %CAN GET THE RAW TIMES FOUND WHERE YOU ALSO CAN USE THE
            %FUNCTION. DONT USE THIS FUNCTION UNLESS IN REFERENCE CLASS
            this.rawtT = this.rawtT/this.reference.nueeRef*nueeRefNew; 
            this.rawtn = this.rawtn/this.reference.nueeRef*nueeRefNew; 
            this.rawtZ = this.rawtZ/this.reference.nueeRef*nueeRefNew; 
            this.rawtE = this.rawtE/this.reference.nueeRef*nueeRefNew;
            if ~isempty(this.rawspecies)
                this.rawspecies.times = this.rawspecies.times/this.reference.nueeRef*nueeRefNew;
            end
        end
    end
    
end
