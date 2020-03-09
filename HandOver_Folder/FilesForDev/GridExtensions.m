    %Check for grid extension
    if ss.useAutomaticGridExtension && isStepToCheckGrid(iteration)

        [doExtend,usefulFrac] = CheckForGridExtension(fMinus1,ss.Ny,ss.Nxi,ss.usefulThreshold);

        if doExtend
            NyOld = ss.Ny;
            NxiOld = ss.Nxi;
            [fMinus1,y,ddy,d2dy2,yWeights,ss,matrixSize,fMinus2] = ...
                DoExtendGrid(fMinus1,ss,gridNormGridPoint,gridNormVal,...
                sp.maxYForGridGrowth,fMinus2);

            if ss.updateSettingsObject
                settings.Nxi = ss.Nxi;
                settings.Ny = ss.Ny;
                settings.yMax = ss.yMax;
                settings.yMaxForPlot = ss.yMax;
                settings.useAutomaticGridExtension = ss.useAutomaticGridExtension;
            end

            %Re-initialize quantities related to the grid and
            %the avalanche source
            nr_mask = (deltaRef*y) > nr_threshhold_w; %%%A: non-runaway removing mask

            switch ss.runawayRegionMode
                case 0 %Isotropic
                case 1 %xi-dependent y_c
                    densityWeightsMat = spdiags((y2.*yWeights)',0,ss.Ny,ss.Ny);
                    currentWeightsMat = spdiags((xy2.*yWeights)',0,ss.Ny,ss.Ny);
                    energyWeightsMat = spdiags((energyMoment.*yWeights)',0,ss.Ny,ss.Ny);
                    [cumintMat,xiGrid] = GenerateCumIntMat(ss.nXi,ss.nPointsXiInt);
                otherwise
                    error('Invalid runaway region mode');
            end

            if idTail %Update the definition of "fast" particle
                tail_mask = y > y(idTail);
            else
                tail_mask = zeros(size(y));
            end

            switch this.eqSettings.sourceMode
                case 3
                    [cHSourcePreFactor,idYMaxForCHSource,cHSourceMatrix,...
                        cHMomentumVector,sourceVector,preFactorConstant,...
                        sigmaIntVector,sinkIndices,gammaM] = ...
                        PreliminariesForChiuHarveySource(this.eqSettings.sourceMode,y,y2,...
                        x,ss.Nxi,deltaRef,sp.nueeRef,ss.nRef,ss.yCutSource);
                    idYFirstAboveCrit = find(mask,1);
                    wMinCH = 2 * sqrt(gammaC*(gammaC-1));
                    intMask = w >= wMinCH;
                    idYMinCHInt = find(intMask,1);
                    cHIntVector = this.plasma.neTotalOverneFree(iteration)*yWeights.*cHMomentumVector.*intMask;
                otherwise
                    %Nothing to do
            end

            if ss.keepTrackOfExternalRunaways
                errorFuncs = erf(x/sp.veBars(1));
                %With scalar parameters its enough to do this once, rather than inside
                %GetFlowThroughGridBoundary, which saves time (erf is quite slow).

                yBoundaryForExternal = round(0.95*y(end));
            end


            %We also need to extend the matrix where we save the
            %output steps, and all the saved distributions, to
            %match the new grid
            fWholeNew = zeros(matrixSize,numel(indicesToReturn));
            fNew = zeros(ss.Ny,nTimeSteps);
            addNy = ss.Ny-NyOld;
            addNxi = ss.Nxi-NxiOld;
            for iExtend = 1:iSaveIndex
                fWholeNew(:,iExtend) = ...
                    ExtendDist(fWhole(:,iExtend),NxiOld,NyOld,addNxi,addNy);
            end
            fNew(1:NyOld,1:iteration) = f(1:NyOld,1:iteration);
            f = fNew;
            fWhole = fWholeNew;
            clear fNew fWholeNew


            ss.yMaxForPlot = ss.yMax;
            f0 = fMinus1(1:ss.Ny);
            gridWasExtended = true;
            gridExtensions(iteration) = 1;
        end
    end
