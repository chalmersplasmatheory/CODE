classdef EfieldOperator < ImplicitOperator
    %EFIELDOPERATOR Electric field operator. It is a tridiagonal.

    properties (SetAccess = private)
        RelLimitUsed = 0
        yUsed = 0
        EHatUsed = 0
        NxiUsed = 0
    end
    
    methods
        function this = EfieldOperator(state,eqSettings)
            %EFIELDOPERATOR constructor. Calls parent to set grid, plasma
            %and eqSettings in this object. For derivation, see
            this@ImplicitOperator(state,eqSettings);
        end

        function matrixHasChanged = generateOperatorMatrix(this,iteration)
            % set values            
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
                y = this.state.momentumGrid.x;
                ddy = this.state.momentumGrid.ddx;
            else
                y = this.state.momentumGrid.y;
                ddy = this.state.momentumGrid.ddy;
            end
            EHat = this.state.physicalParams.EHats(iteration);
            Ny = this.state.momentumGrid.Ny;
            Nxi = this.state.momentumGrid.Nxi;

            %Check if matrix needs to be built
            if this.RelLimitUsed == takeNonRelLimit && ...
                    isequal(this.yUsed,y) && ...
                    this.EHatUsed == EHat && ...
                    this.NxiUsed == Nxi
                matrixHasChanged = 0;
                return
            end
            
            disp('Efield Builds');
            %Build matrix
            this.predictedNNZ = 2 * Nxi * Ny + 2 * Nxi * nnz(abs(ddy));
            this.resetSparseCreator();
            
            %TODO Row range should vary depending on boundary condition?
            for L=0:(Nxi-1)
                rowIndices = L*Ny + (2:(Ny-1));

                % Add electric field term,
                % Sub-diagonal term in L:
                ell = L-1;
                if ell >=0
                    columnIndices = ell*Ny + (1:Ny);
                    littleMatrix = EHat*L/(2*L-1)*ddy - spdiags((EHat*(L-1)*L/(2*L-1)./y)',0,Ny,Ny);
                    this.addSparseBlock(rowIndices, columnIndices, littleMatrix(2:(Ny-1), 1:Ny))
                end

                % Add electric field term,
                % Super-diagonal term in L:
                ell = L+1;
                if ell<Nxi
                    columnIndices = ell*Ny + (1:Ny);
                    littleMatrix = EHat*(L+1)/(2*L+3)*ddy + spdiags((EHat*(L+1)*(L+2)/(2*L+3)./y)',0,Ny,Ny);
                    this.addSparseBlock(rowIndices, columnIndices, littleMatrix(2:(Ny-1), 1:Ny))
                end
            end

            this.operatorMatrix = -this.createSparse();
            this.clearSparseResidue();
            
            %save values used
            this.RelLimitUsed = takeNonRelLimit;
            this.yUsed = y;
            this.EHatUsed = EHat;
            this.NxiUsed = Nxi;
            
            matrixHasChanged = 1;
        end
        
    end
end

