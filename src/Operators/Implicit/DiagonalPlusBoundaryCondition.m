classdef DiagonalPlusBoundaryCondition < ImplicitOperator
    %DIAGONALPLUSBOUNDARYCONDITION builds the diagonal and boundary condition matrix
    %in CODE version 1.03 CODEtimeDependent(). Basically does the BC and
    %diag matrix under PreliminariesForMatrix
    %done under the switch variable CollisionalOperator, generate matrix is
    %the same as BuildMatrix function

    % Properties used for the given operator
    properties
        NxiUsed = 0
        ddyUsed = 0
        RelLimitUsed = 0
        yMaxBC = 0
    end
    
    methods
        function this = DiagonalPlusBoundaryCondition(state,eqSettings)
            %DIAGONALPLUSBOUNDARYCONDITION Construct an instance of this class
            this@ImplicitOperator(state,eqSettings);
        end

        function matrixHasChanged = generateOperatorMatrix(this,iteration)
            %generateOperatorMatrix generates the diagonalplusboundarycondition matrix
            %%% Build a matrix with the boundary conditions and a diagonal %%

            Ny = this.state.momentumGrid.Ny;
            Nxi = this.state.momentumGrid.Nxi;
            yMax = this.state.momentumGrid.yMax;
            matrixSize = this.state.momentumGrid.matrixSize;

            switch this.eqSettings.collisionOperator
                case {0,2,3,4}
                    if this.eqSettings.useNonRelativisticCollOp
                        takeNonRelLimit = true;
                    else
                        takeNonRelLimit = false;
                    end
                case {1}
                    takeNonRelLimit = true;
                otherwise
                    error('You have chosen an invalid collision operator.');
            end

            if takeNonRelLimit
                ddy = this.state.momentumGrid.ddx;
            else
                ddy = this.state.momentumGrid.ddy;
            end

            if this.NxiUsed == Nxi && ...
                    isequal(this.ddyUsed, ddy) && ...
                    this.RelLimitUsed == takeNonRelLimit && ...
                    this.yMaxBC == this.state.momentumGrid.yMaxBoundaryCondition
                matrixHasChanged = 0;
                return
            end

            disp('Diag&BC builds');
            %Build Matrix
            
            % Predict roughly how many nonzero elements will be in the matrix of
            % boundary conditions. This speeds up the code by eliminating the need
            % to reallocate memory during matrix construction:
            this.predictedNNZ = 2*(Nxi + matrixSize) + Ny;
            this.resetSparseCreator()

            indices = 2:matrixSize;
            this.addToSparse(indices, indices, ones(size(indices)))

            % For the special point at y=0, apply Neumann condition dF/dy=0 for L=0:
            this.addSparseBlock(1,1:Ny,ddy(1,:));

            for L = 1:(Nxi-1)
                % Add boundary conditions at yMax:
                index = L*Ny + Ny;
                this.addToSparse(index, index, 1);
                %Add boundary conditions at y=0:
                index = L*Ny + 1;
                this.addToSparse(index, index, 1);
            end

            % Impose boundary conditions at y=yMax:
            if this.state.momentumGrid.yMaxBoundaryCondition==2
                % Add Robin boundary condition dF/dy + (2/y)*F = 0 at yMax:
                for L=0:(Nxi-1)
                    rowIndex = L*Ny + Ny;
                    columnIndices = L*Ny + (1:Ny);
                    boundaryCondition = ddy(Ny,:);
                    boundaryCondition(Ny) = boundaryCondition(Ny) + 2/yMax - 1; %Subtract 1 since we already put a 1 on the diagonal.
                    this.addSparseBlock(rowIndex, columnIndices, boundaryCondition)
                end
            end

            this.operatorMatrix = this.createSparse();
            this.clearSparseResidue();
            
            %save used values for given matrix
            this.NxiUsed = Nxi;
            this.ddyUsed = ddy;
            this.RelLimitUsed = takeNonRelLimit;
            this.yMaxBC = this.state.momentumGrid.yMaxBoundaryCondition;
            matrixHasChanged = 1;
        end
    end
end

