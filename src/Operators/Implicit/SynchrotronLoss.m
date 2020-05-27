classdef SynchrotronLoss < ImplicitOperator
    %SynchrotronLoss
    
    % Properties used for the given operator
    properties (SetAccess = private)
        nueeBarUsed = 0
        BHatUsed = 0
        NxiUsed = 0
        yUsed = 0
        deltaRefUsed = 0
        gammaUsed = 0
        yMaxBCUsed = 0
        useFullSynchOpUsed = 0;
    end
    
    methods
        function this = SynchrotronLoss(state,eqSettings)
            %SynchrotronLoss object to calculate syncrotron loss term in
            %code run
            this@ImplicitOperator(state,eqSettings)
        end
    end

    methods

        function matrixHasChanged = generateOperatorMatrix(this,iteration)
            %Initialize values
            
            nueeBar = this.state.physicalParams.nueeBars(iteration);
            BHat = this.state.physicalParams.BHatRef;
            Ny = this.state.momentumGrid.Ny;
            Nxi = this.state.momentumGrid.Nxi;
            matrixSize = Nxi*Ny;
            
            y = this.state.momentumGrid.y;
            y2 = this.state.momentumGrid.y2;
            ddy = this.state.momentumGrid.ddy;
            deltaRef = this.state.reference.deltaRef;
            gamma = this.state.momentumGrid.gamma;
            x = this.state.momentumGrid.x;
            %
            if this.nueeBarUsed == nueeBar && ...
                    isequal(this.BHatUsed, BHat) && ... % Ugly hack, sometimes BHat is empty, so [] == [] is the same as [], which is treated as false.
                    this.NxiUsed == Nxi && ...
                    isequal(this.yUsed,y) && ...
                    this.deltaRefUsed == deltaRef && ...
                    isequal(this.gammaUsed, gamma) && ...
                    this.yMaxBCUsed == this.state.momentumGrid.yMaxBoundaryCondition &&...
                    this.useFullSynchOpUsed == this.eqSettings.useFullSynchOp
                disp('not here')
                    matrixHasChanged = 0;
                return
            end
            if ~isempty(BHat)
                disp('Synch builds');
                %Matrix Building
                this.predictedNNZ = 3 * Nxi * (Ny + nnz(abs(ddy)));
                this.resetSparseCreator()

                if this.state.momentumGrid.yMaxBoundaryCondition == 3
                    rowRange = 2:Ny;
                else
                    rowRange = 2:(Ny-1);
                end

                %%% This is the synchrotron radiation loss term %%%     
                rootFactor = x./y; % 1/sqrt(1+(deltaRef*y).^2);
                preFactor = BHat^2*rootFactor; %The sign is reversed compared to the derivation!                
                if this.eqSettings.useFullSynchOp
                    bracket = diag(4*deltaRef^2*preFactor.*y2) + diag(preFactor.*y.*(1+deltaRef^2.*y2))*ddy;
                else
                    bracket = diag(preFactor.*y.*(1+deltaRef^2.*y2))*ddy;
                end
                pitchADiag = @(L) 1 - (L+1)^2/((2*L+3)*(2*L+1)) - L*L/(4*L*L-1);
                pitchBDiag = @(L) ( L*L/(2*L-1) - (L*(L+1)^2)/((2*L+1)*(2*L+3)) - (L^3)/((2*L-1)*(2*L+1)) );% = L*L*(L+1)/(4*L*L-1) - L*(L+1)^2/((2*L+1)*(2*L+3));
                pitchASub = @(L) -L*(L-1)/((2*L-1)*(2*L-3));
                pitchBSub = @(L) -L*(L-1)*(L-2)/((2*L-1)*(2*L-3));
                pitchASup = @(L)  -(L+1)*(L+2)/((2*L+3)*(2*L+5));
                pitchBSup = @(L)  (L+1)*(L+2)/(2*L+3) - (L+1)*(L+2)^2/((2*L+3)*(2*L+5)); % = (L+1)*(L+2)*(L+3)/((2*L+3)*(2*L+5));

                for L = 0:(Nxi-1)
                   %%%Diagonal term
                   if this.eqSettings.useFullSynchOp
                       smallMatrix = diag(2*preFactor); 
                   else
                       smallMatrix = [];
                   end
                   smallMatrix =  smallMatrix + bracket*pitchADiag(L) - diag(preFactor*pitchBDiag(L));      
                   rowIndices = L*Ny + rowRange;
                   columnIndices = L*Ny + (1:Ny);
                   this.addSparseBlock(rowIndices, columnIndices, smallMatrix(rowRange, 1:Ny))       
                   %%%Sub-diagonal term
                   ell = L-2;
                   if ell >= 0                
                       smallMatrix = bracket*pitchASub(L) - diag(preFactor*pitchBSub(L));     
                       columnIndices = ell*Ny + (1:Ny);
                       this.addSparseBlock(rowIndices, columnIndices, smallMatrix(rowRange, 1:Ny))
                       if ell==0
                           this.addSparseBlock(rowIndices, matrixSize, smallMatrix(1+rowRange, 1))
                       end
                   end
                   %%%Super-diagonal term
                   ell = L+2;
                   if ell < Nxi                
                       smallMatrix = bracket*pitchASup(L) - diag(preFactor*pitchBSup(L));
                       columnIndices = ell*Ny + (1:Ny);            
                       this.addSparseBlock(rowIndices, columnIndices, smallMatrix(rowRange, 1:Ny))
                   end
                end


                this.operatorMatrix = this.createSparse();
                this.clearSparseResidue()

                this.operatorMatrix = -(1/nueeBar)*this.operatorMatrix;

                %save values relevant for knowing whether to rebuild matrix or
                %not
                this.nueeBarUsed = nueeBar;
                this.BHatUsed = BHat;
                this.NxiUsed = Nxi;
                this.yUsed = y;
                this.deltaRefUsed = deltaRef;
                this.gammaUsed = gamma;
                this.yMaxBCUsed = this.state.momentumGrid.yMaxBoundaryCondition;
                this.useFullSynchOpUsed = this.eqSettings.useFullSynchOp;
                
                matrixHasChanged = 1;
            else
                this.operatorMatrix = sparse(matrixSize,matrixSize);
                matrixHasChanged = 0;
            end
        end

    end
end

