classdef Solver < handle
    % Abstract solver class structuring a solver class to be used in CODE

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Objects describing plasma:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   state - State object containing physical information about the
    %            plasma such as temperature, Efield, Bfield and about what momentum grid
    %          is used. This is a handle class shared by all operator
    %          objects
    %   implicitOp - cell containing all implicit operators used in the run
    %              all objects in cell should extend ImplicitOperator class
    %   explicitOp - cell containing all explicit operators used in the run
    %              all objects in cell should extend ExplicitOperator class
    %   sources - cell containing all sources for the run, care these are
    %               not decided how they should be structured inside
    %               the code yet
    %   eqSettings - EquationSettings containing what settings for the
    %               equation to use
    properties (SetAccess = protected)
        % constant vars under the run
        state
        eqSettings
        implicitOperators = {}
        explicitOperators = {}
        sources = {}
    end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%                     %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%     Constructor     %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%                     %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%                     %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function this = Solver(state, eqSettings)
            if ~isa(eqSettings,'EquationSettings')
                error('Second argument should be a EquationSettings object.')
            end
            this.eqSettings = eqSettings;
            this.updateOperators;
            if isa(state,'State') 
                this.state = state;
            else
                error('The constructor only accepts State object.');
            end
        end

    end


       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%     update          %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%     operators       %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%                     %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function updateOperators(this)
            ik = 1; sk = 1; ek = 1;
            for k = 1:length(this.eqSettings.operators)
                if isa(this.eqSettings.operators{k},'ImplicitOperator')
                    this.implicitOperators{ik} = this.eqSettings.operators{k};
                    ik = ik+1;
                elseif isa(this.eqSettings.operators{k},'ExplicitOperator')
                    this.explicitOperators{ek} = this.eqSettings.operators{k};
                    ek = ek+1;
                elseif isa(this.eqSettings.operators{k},'Source')
                    this.sources{sk} = this.eqSettings.operators{k};
                    sk = sk+1;
                elseif isempty(this.eqSettings.operators{k})
                else
                    error('Unknown operator type in passed eqSettings.')
                end
            end
        end
    end
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%     time            %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%     stepping        %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%     Functions       %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%                     %%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods (Abstract)
            takeTimeSteps(this)
   end

end

