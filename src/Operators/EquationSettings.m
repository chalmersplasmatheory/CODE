classdef EquationSettings <handle
    %EQUATIONSETTINGS settings class defining what numerical and analytical
    %approximations and/or features to use.

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Collision operators and conservational properties:   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % collisionOperator - specifies the type of collision operator to use
    %       0 = Fully relativistic test-particle collision operator.
    %           NOT momentum-conserving. Relativistic Fokker-Planck
    %           coefficients from OJ Pike & SJ Rose, Phys Rev E 89 (2014).
    %           (Default)
    %       1 = Non-relativistic, linearized, momentum-conserving. Includes
    %           modified Rosenbluth potentials.
    %       2 = Use the non-relativistic field particle term from 1, together
    %           with the relativistic test particle term from 4.
    %           Generally works better than 1, as the tail quickly becomes
    %           relativistic. The field-particle term only affects the
    %           bulk, so the non-relativistic approximation is usually good
    %           there.
    %       3 = Use an ad-hoc relativistic generalization of the
    %           field-particle term from 1, together with the test-particle
    %           terms from 4. Generally appears to work better than 1 & 2,
    %           although the physics basis is shaky. !!!Experimental!!!
    %       4 = Relativistic collision operator for non-relativistic
    %           thermal speeds. Agrees with 0 when T << 10 keV.
    %           See the CODE-paper for details.
    % useNonRelativisticCollOp - use the non-relativistic limit of the
    %                            collision operator, for instance for
    %                            benchmarking hot-tail generation against
    %                            Smith et al. (2005), or Smith and
    %                            Verwichte (2008). (Default: false)

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Avalanche (secondary runaway) source term:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sourceMode - the type of avalanche source to use
    %       0 = do not include a source (Default)
    %       1 = include a Rosenbluth-Putvinskii-like source
    %       2 = the same sorce, but using the number of "fast" particles,
    %           rather than the number of runaways, to determine the
    %           source magnitude. Useful in runaway decay and/or when E<E_c,
    %           (when there are no runaways by definition).
    %           DISCLAIMER: This is not necessarily consistent, and is
    %           sensitive to the definition of a fast particle!
    %       3 = include a Chiu-Harvey-like source
    %       4 = include a Chiu-Harvey-like source which takes the
    %           pitch-dependence of the distribution into account (derived
    %           by Ola). Should be more correct than 3.
    %       5 = The same as 4, but including a sink term to conserve
    %           particle energy and momentum.
    % fastParticleDefinition - relevant when sourceMode = 2
    %       0 = Do not calculate the fast particle content, is not
    %           compatible with sourceMode 2
    %       1 = Based on where the "tail" "begins" (where f/f_M > some threshold)
    %			Note that the shape of the source in momentum space is only
    %			recalculated if the electric field changes, and may at times not
    %			be consistent with the beginning of the tail.
    %       2 = Based on relative speed (v/v_e > some threshold)
    %       3 = Based on absolute speed (v/c > some threshold)
    %       4 = Selects the most generous of cases 2,3 and the critical speed for
    %           runaways. Useful when you want to follow n_r closely, but still
    %           have fast particle when E is close to (or below) E_c.
    %            (Default)
    % tailThreshold - for fastParticleDefinition = 1 (Default: 1e3)
    % relativeSpeedThreshold - for fastParticleDefinition = 2 (Default: 10)
    % absoluteSpeedThreshold - for fastParticleDefinition = 3 (Default:
    % 0.25)
    % yCutSource - Arbitrary cutoff needed for source 5. Set to empty to use
    %              yCutSource = y_crit (currently only works with constant
    %              plasma parameters) Can also be used for source 3 and 4
    %              (cf. sourceMode 1 vs 2) (Default: 5)
    % runawayRegionMode - determines how to define the runaway region and
    %                     how to calcualate moments of f in that region.
    %			Is interpreted as 0 no mather what in sourceMode 2.
    %       0 = Isotropic (y_c=y_c(xi=1) for all xi) (Default)
    %       1 = xi-dependent y_c -- more correct. Most relevant for
    %           hot-tail scenarios
    % nPointsXiInt - resolution parameter for runawayRegionMode = 1
    %                (Default: 200)


    % %%%%%%%%%%%%%%%%
    % Miscellaneous:
    % %%%%%%%%%%%%%%%%

    % enforceDensityConservation - enables an adaptive particle source term
    %                              which keeps the density constant by
    %                              adding/removing partiles to/from the
    %                              bulk. Useful in runs with many time
    %                              steps, when the step-to-step numerical
    %                              losses accumulate to a large error.
    %                               (Default: true)
    % useScreening - take screening into account
    %       0 = Ignore screening (complete screening, correct at low
    %       energy) (Default)
    %       1 = Fokker--Planck collision operator, TF model (works for all species)
	%		2 = Fokker--Planck collision operator, TF-DFT model (works for some
	%			ionization degrees of argon and Xe)
	%		3 = Fokker--Planck collision operator, full DFT model (only works for Ar+ so far)
	% 		4 = Boltzmann collision operator, TF model
	%		5 = Boltzmann collision operator, TF-DFT model
	%		6 = Boltzmann collision operator, full DFT model
    %       7 = Full penetration (no screening, Z0^2-> Z^2, correct in the
    %       limit p-> infinity)
    %
    % useInelastic: take inelastic collisions with bound electrons into
    %       account
    %       0 = Ignore inelastic collisions (Default)
    %       1 = Bethe formula
    %       2 = RP
    %
    % useEnergyDependentLnLambdaScreening: logarithmic enhancement of the
    %                                       Coulomb logarithm
    %       0 = No (standard formula) (Default)
    %       1 = Yes (change in the collision operator). Interpolated Landau form.
    %       2 = Enhance by an ad-hoc factor of 1.3.
    %
    % NyInterp - number of grid points in the interpolated grid necessary
    %            for the Boltzmann-based bremsstrahlung and knock-on
    %            operators (bremsMode>=2, sourceMode = 5)
    % bremsMode - include an operator for bremsstrahlung radiation reaction
    %       2 = Boltzmann op., full (both low and high k, angular-dependent
    %           cross section)
    %       3 = Boltzmann op., low and high k, no angular dependence
    %       4 = Boltzmann op., only high k, no angular dependence
    % useFullSynchOp - determines whether to use the full, correct,
    %                  synchrotron radiation reaction operator, or just
    %                  parts of it (for benchmarking to Andersson, et al., PoP, 2001)
    properties (SetAccess = protected)
        % Collision Operator
        collisionOperator = 0;
        useNonRelativisticCollOp = false;

        %Avalanche source
        sourceMode = 0;
        fastParticleDefinition = 4;
        tailThreshold = 1e3;
        relativeSpeedThreshold = 10;
        absoluteSpeedThreshold = 0.25;
        yCutSource = 5
        runawayRegionMode = 0;
        nPointsXiInt = 200;

        %Misc
        enforceDensityConservation = true
        useScreening = 0;
        useInelastic = 0;
        useEnergyDependentLnLambdaScreening = 0;
        bremsMode = 0;
        useFullSynchOp = 1; 
        NyInterp = 1000;
        %Operators and state
        operators = cell(6,1);
        state
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                %%%%%%%%%%%%%%%
    %%%%%%%%%%%%    Version     %%%%%%%%%%%%%%%
    %%%%%%%%%%%%                %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
        VERSION = 1.0
    end

    methods
        function this = EquationSettings(state,varargin)
            %EQUATIONSETTINGS Construct an instance of this class
            %To initiate properties, write "property name", property
            if (isa(state,'State') && state.VERSION == this.VERSION)
                this.state = state;
            end
            for k = 2:2:nargin
                if strcmp('operators', varargin{k})
                   error('Setting the operators by the user is not supported yet.')
                end
                this.(varargin{k}) = varargin{k+1};
            end
            this.setOperators();
        end
        
        function setOperators(this)
            if isempty(this.operators{1})
                this.operators{1} = DiagonalPlusBoundaryCondition(this.state,this);
                this.operators{2} = CollisionOperator(this.state,this);
            end
            
            if isempty(this.operators{3})
                this.operators{3} = EfieldOperator(this.state,this);
            end
            
            if isempty(this.operators{4})
                this.operators{4} = SynchrotronLoss(this.state,this);
            end

            switch this.sourceMode
                case 0 
                    this.operators{5} = {};
                case {1,2} 
                    if ~isa(this.operators{5}, 'RosenbluthSource')
                        this.operators{5} = RosenbluthSource(this.state,this);
                    end
                case {3} 
                    if ~isa(this.operators{5}, 'ChiuHarveySource')
                        this.operators{5} = ChiuHarveySource(this.state,this);
                    end
                case {4,5}
                   if ~isa(this.operators{5}, 'ImprovedChiuHarveySource')
                      this.operators{5} = ImprovedChiuHarveySource(this.state,this); 
                  end
              otherwise
                  error('Unknown sourceMode')
            end
            switch this.bremsMode
                case {0,1}
                    this.operators{6} = {};
                case {2,3,4}
                    this.operators{6} = Bremsstrahlung(this.state,this);
                otherwise
                    error('Unknown bremsMode')
            end
            
        end

    end

    methods 
        
        function setcollisionOperator(this,value)  
 		    this.collisionOperator = value;
 		end

        function setuseNonRelativisticCollOp(this,value)  
 		    this.useNonRelativisticCollOp = value;
 		end

        function setsourceMode(this,value)  
 			this.sourceMode = value;
            this.setOperators();
 		end

        function setfastParticleDefinition(this,value)  
 			this.fastParticleDefinition = value;
 		end

        function settailThreshold(this,value)  
 			this.tailThreshold = value;
 		end

        function setrelativeSpeedThreshold(this,value)  
 			this.relativeSpeedThreshold = value;
 		end

        function setabsoluteSpeedThreshold(this,value)  
 			this.absoluteSpeedThreshold = value;
 		end

        function setyCutSource(this,value)  
 			this.yCutSource = value;
 		end

        function setrunawayRegionMode(this,value)  
 			this.runawayRegionMode = value;
 		end

        function setnPointsXiInt(this,value)  
 			this.nPointsXiInt = value;
 		end

        function setenforceDensityConservation(this,value)  
 			this.enforceDensityConservation = value;
 		end

        function setuseScreening(this,value)  
 			this.useScreening = value;
 		end

        function setuseInelastic(this,value)  
 			this.useInelastic = value;
 		end

        function setuseEnergyDependentLnLambdaScreening(this,value)  
 			this.useEnergyDependentLnLambdaScreening = value;
 		end

        function setbremsMode(this,value)  
 			this.bremsMode = value;
            this.setOperators();
 		end
        
        function setNyInterp(this,value)
            this.NyInterp = value;
        end

        function setuseFullSynchOp(this,value)
            this.useFullSynchOp = value;
        end

        function mimic(this,eqSettings)
            %copies all settings from another EquationSettings object
            if ~isa(eqSettings,'EquationSettings')
                error('Can only mimic settings from another settings object');
            end
            if eqSettings.VERSION ~= this.VERSION
                error('Cannot mimic EquationSettings of different versions')
            end
			this.collisionOperator = eqSettings.collisionOperator; 
			this.useNonRelativisticCollOp = eqSettings.useNonRelativisticCollOp; 
			this.sourceMode = eqSettings.sourceMode; 
			this.fastParticleDefinition = eqSettings.fastParticleDefinition; 
			this.tailThreshold = eqSettings.tailThreshold; 
			this.relativeSpeedThreshold = eqSettings.relativeSpeedThreshold; 
			this.absoluteSpeedThreshold = eqSettings.absoluteSpeedThreshold; 
			this.yCutSource = eqSettings.yCutSource; 
			this.runawayRegionMode = eqSettings.runawayRegionMode; 
			this.nPointsXiInt = eqSettings.nPointsXiInt; 
			this.enforceDensityConservation = eqSettings.enforceDensityConservation; 
			this.useScreening = eqSettings.useScreening; 
			this.useInelastic = eqSettings.useInelastic; 
			this.useEnergyDependentLnLambdaScreening = eqSettings.useEnergyDependentLnLambdaScreening; 
			this.bremsMode = eqSettings.bremsMode;
            if this.state == eqSettings.state
                for i = 1:length(eqSettings.operators)
                    if ~isempty(eqSettings.operators{i})
                        this.operators{i} = copy(eqSettings.operators{i});
                    end
                end
            else
                this.setOperators()
            end
        end
    end
end

