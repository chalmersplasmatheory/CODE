classdef SolverSteadyState
    %SOLVERSTEADYSTATE solves the SteadyState equation on the form: MF = S
    %where the source S is on the form S = a exp(-y^2)

    properties
        state
        eqSettings
        operator
    end

    methods
        function this = SolverSteadyState(state, eqSettings)
            %SOLVERSTEADYSTATE Construct an instance of this class
            %   Detailed explanation goes here
            this.state = state;
            this.eqSettings = eqSettings;
        end

        function [F, a] = findSteadyState(this)
            % calculates the steady state solution to operator f = a * exp(-y2)
            % and interprets the a constant as the rate of change of
            % runaways
            this.grid.interpolatePhysicalParams();
            this.plasma.timeNormalize();
            this.plasma.calcDepParams();
            this.grid.initializeGrid();
            
            
            this.operator.generateOperatorMatrix(1);
            sl = SynchrotronLoss(this.grid,this.plasma,this.eqSettings);
            sl.generateOperatorMatrix(1);
            [factor_L, factor_U, factor_P, factor_Q] = lu(-this.operator.operatorMatrix-sl.operatorMatrix);
            rhs = zeros(this.grid.matrixSize,1);
            rhs(1:this.grid.Ny) = exp(-this.grid.y2);
            rhs(this.grid.Ny) = 0;
            Fda = factor_Q * (factor_U \ (factor_L \ (factor_P * rhs)));
            a = Fda(1);
            F = Fda/a;
        end
        
        
        function [F, a] = steadyState2(this)

            % Init operators
            d = DiagonalPlusBoundaryCondition(this.state,this.eqSettings);
            d.generateOperatorMatrix(1);
            co = CollisionOperator(this.state,this.eqSettings);
            co.generateOperatorMatrix(1);
            s = SynchrotronLoss(this.state,this.eqSettings);
            s.generateOperatorMatrix(1);
            ef = EfieldOperator(this.state,this.eqSettings);
            ef.generateOperatorMatrix(1);
            


            M = (s.operatorMatrix+co.operatorMatrix+ef.operatorMatrix);

            rhs = zeros(this.state.momentumGrid.Ny * this.state.momentumGrid.Nxi, 1);
            rhs(1:this.state.momentumGrid.Ny) = exp ( -this.state.momentumGrid.y .^ 2 );

            % Dirichlet condition.
            rhs(1) = 0;
            M(1, 1:5) = -this.state.momentumGrid.ddy(1, 1:5);

            % Neumann at y = 0
            for L = 1:(this.state.momentumGrid.Nxi - 1)
                M(L * this.state.momentumGrid.Ny + 1, L * this.state.momentumGrid.Ny + 1) = 1;
            end

            % Neumann at y = y_max (I don't really like this).
            for L = 0:(this.state.momentumGrid.Nxi - 1)
                M(L * this.state.momentumGrid.Ny + this.state.momentumGrid.Ny, L * this.state.momentumGrid.Ny + this.state.momentumGrid.Ny) = 1;
            end

            result = -(M \ rhs);

            if this.eqSettings.collisionOperator == 2
                %minus sign from that there is a minus in front of -M \rhs
                result(1:this.state.momentumGrid.Ny) = result(1:this.state.momentumGrid.Ny) - addHeat_ParticleSources(result, this.state.momentumGrid.Ny,this.state.momentumGrid.y2,this.state.momentumGrid.gamma,this.state.momentumGrid.yWeights, this.state.physicalParams.nBars(1), this.state.physicalParams.veBars3(1), this.state.physicalParams.veBars2(1), this.state.physicalParams.deltaRef^2);
            end
            F = result / result(1);
            a = -1 / result(1);
        end
        
    end
end

