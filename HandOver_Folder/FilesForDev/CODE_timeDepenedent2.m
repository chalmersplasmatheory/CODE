function varargout = CODE_timeDependent(varargin)

if sp.maxYForGridGrowth <= this.plasma.yMax
    if this.plasma.useAutomaticGridExtension %Print info in case we thought we could use grid extension
        fprintf('\nNo further grid extensions are allowed by the value of the externalBoundaryMeV property\n');
    end
    this.plasma.useAutomaticGridExtension = false;
    if this.plasma.updateSettingsObject
        settings.useAutomaticGridExtension = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the building of the matrix and initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize arrays for building the sparse matrix:

    nr_mask = (deltaRef*y) > nr_threshhold_w; %%%A: non-runaway removing mask

switch ss.runawayRegionMode
    case 0 %Isotropic
    case 1 %xi-dependent y_c
        densityWeightsMat = spdiags((y2.*yWeights)',0,ss.Ny,ss.Ny);
        currentWeightsMat = spdiags((xy2.*yWeights)',0,ss.Ny,ss.Ny);
        energyWeightsMat = spdiags((energyMoment.*yWeights)',0,ss.Ny,ss.Ny);
        [cumintMat,xiGrid] = GenerateCumIntMat(ss.Nxi,ss.nPointsXiInt);
    otherwise
        error('Invalid runaway region mode');
end
switch this.eqSettings.sourceMode
    case 3
        [cHSourcePreFactor,idYMaxForCHSource,cHSourceMatrix,cHMomentumVector,...
            sourceVector,preFactorConstant,sigmaIntVector,sinkIndices,gammaM] = ...
            PreliminariesForChiuHarveySource(this.eqSettings.sourceMode,y,y2,x,ss.Nxi,deltaRef,...
            sp.nueeRef,ss.nRef,ss.yCutSource);
    otherwise
        %Nothing to do
end
if ss.keepTrackOfExternalRunaways

    yBoundaryForExternal = round(0.95*y(end));

    errorFuncs = erf(x/sp.veBars(1));
    %With scalar parameters its enough to do this once, rather than inside
    %GetFlowThroughGridBoundary, which saves time (erf is quite slow).
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of time-advance loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diffOut = zeros(size(y));

%Calculate the boundary for the fast-particle definition which we need to
%build the avalanche source in the first time step
if this.eqSettings.sourceMode == 2
    [idTail,tail_mask] = DetermineFastParticleRegion(ss,x,y,fMinus1,...
        nBar,veBar3,deltaRef,0);
end

for iteration = 1:nTimeSteps

    %%%%% Compute moments of the distribution in the previous time step %%%
    % Save also the f0 part of every time step
    f(:,iteration) = soln(1:ss.Ny);

    %Update runaway region
    if sp.parametersHaveChanged(iteration) || gridWasExtended
        %Update the runaway region and source (nr_threshhold_w &
        %wMinForSource depend on E/E_c)
        if sp.noTimeDep && sp.EOverEc(1) > 1   %@@@This should be improved!
            gammaC = sqrt( sp.EOverEc(1)/(sp.EOverEc(1)-1) );
            wMinForSource = sqrt(4*gammaC*(gammaC-1)); % only count runaways with
            % kinetic energy twice the
            % critical runaway energy
            %            wMinForSource = 1/sqrt(sp.EOverEc(1)-1);
        elseif ~sp.noTimeDep && sp.EOverEc(iteration) > 1
            gammaC = sqrt( sp.EOverEc(iteration)/(sp.EOverEc(iteration)-1) );
            wMinForSource = sqrt(4*gammaC*(gammaC-1)); % only count runaways with
            % kinetic energy twice the
            % critical runaway energy
            %            wMinForSource = 1/sqrt(sp.EOverEc(iteration)-1);
        else
            wMinForSource = inf;
            gammaC = inf;
        end
        nr_threshhold_w = wMinForSource;
        nr_mask = w > nr_threshhold_w;

        if this.eqSettings.sourceMode == 2
            %Apply the source in the fast-particle region, not the runaway
            %region.
            if ss.fastParticleDefinition == 4 && iteration == 1
                %In the first time step, the runaway definition might still
                %be the most generous
                mask = (w >= wMinForSource) | tail_mask;
            else
                mask = tail_mask;
            end
        else
            mask = w >= wMinForSource;
        end


        if any(mask ~= oldMask) || iteration == 1  || gridWasExtended
            %Rebuild the source if the mask has changed (more/fewer
            %grid points should be included in the source)
            sourceRecomputations(iteration) = true;


            switch this.eqSettings.sourceMode
                case {1,2}

                case 3
                    %Chiu-Harvey-like source
                    idYFirstAboveCrit = find(mask,1);
                    wMinCH = 2 * sqrt(gammaC*(gammaC-1));
                    intMask = w >= wMinCH;
                    idYMinCHInt = find(intMask,1);
                    cHIntVector = this.plasma.neTotalOverneFree(iteration)*yWeights.*cHMomentumVector.*intMask;

                otherwise
                    %Nothing to do
            end
            oldMask = mask;
        end

        gridWasExtended = false;

        switch ss.runawayRegionMode
            case 0 %Isotropic
                %Already done above
            case 1 %xi-dependent y_c
                xiIntMat = MapToXiIntMat(cumintMat,xiGrid,sp.EOverEc(iteration),y2,deltaRef);
            otherwise
        end
    end

    %%% Compute "number of runaways" and current:
    f0 = fMinus1(1:ss.Ny);
    nInst = (yWeights.*y2) * f0;
    switch ss.runawayRegionMode
        case 0 %Isotropic
            nr = (yWeights.*y2.*nr_mask) * f0;
        case 1 %xi-dependent y_c
            nr = CalculateRunawayMomentXiDependent(ss,fMinus1,xiIntMat,densityWeightsMat);
        otherwise
            error('Not implemented yet')
    end
    nr = 4/sqrt(pi) * nr;
    % Note: nr is really the density of runaways normalized by the
    % instantaneous total electron density (n_r/n). Some prefactors
    % have been cancelled, so nInst can not be used directly as a
    % measure of the density
    nrs(iteration) = max(nr,0); %@@@Investigate the implications for plots, etc of using instantaneous n_r/n
    nes(iteration) = nInst;


    %%% Also compute the number of "fast particles"
    %Periodically find out where the tail begins
    if isStepToCheckGrid(iteration) || iteration == 1
        [idTail,tail_mask] = DetermineFastParticleRegion(ss,x,y,fMinus1,...
            nBar,veBar3,deltaRef,nr_mask);
    end
    nFast = 4/sqrt(pi) * (yWeights.*y2.*tail_mask) * f0 ;
    nfs(iteration) = nFast;

    f1 = fMinus1((ss.Ny+1):(2*ss.Ny));
    f2 = fMinus1((2*ss.Ny+1):(3*ss.Ny));

    %Keep track of particles that leave the grid. If the dist it not close
    %to the grid boundary, we don't have to perform the calculation.
    if ss.keepTrackOfExternalRunaways && ~ss.useAutomaticGridExtension
        [~,usefulFrac] = FindLastUsedGridPoint(fMinus1,ss.Ny,ss.Nxi,ss.usefulThreshold);
    end
    if ss.keepTrackOfExternalRunaways && ~ss.useAutomaticGridExtension && (usefulFrac > 0.95 || nExternal > 0)
        if sp.noTimeDep
            it = 1;
        else
            it = iteration;
            %Not very pretty, but saves time to do it outside
            %GetFlowThroughBoundary (we only need to do it once if
            %parameters don't change).
            errorFuncs = erf(x/sp.veBars(it));
        end
        %%% Compute the increase in the number of external runaways
        flowRate = GetFlowThroughGridBoundary(f0,f1,f2,y,y2,x,x2,ddy,...
            deltaRef,nBar,errorFuncs,yBoundaryForExternal,...
            includeSynchrotronLoss,sp,it);
        nExternal = nExternal + flowRate*dt; %flowRate*dt = n_lost/n_e
        nExternals(iteration) = nExternal;
        externalCurrent = -(3*sqrt(pi))/4*nExternal/deltaRef; % nExternal*c/v_e,ref = nExternal/deltaRef

        totalCurrent = (xy2.*yWeights) * f1 + externalCurrent;

        if totalCurrent == 0
            currents(iteration) = 0;
            totRECurrentDens(iteration) = 0;
            totCurrentDens(iteration) = 0;
        else
            switch ss.runawayRegionMode
                case 0
                    runawayCurrent = (xy2.*yWeights.*nr_mask) * f1;
                case 1
                    runawayCurrent = CalculateRunawayMomentXiDependent(ss,fMinus1,xiIntMat,currentWeightsMat);
                otherwise
                    error('Invalid mode');
            end
            fastCurrent = (xy2.*yWeights.*tail_mask) * f1;
            currents(iteration) = runawayCurrent/totalCurrent; % Constants and normalization have been cancelled!
            fastCurrents(iteration) = fastCurrent/totalCurrent; % Constants and normalization have been cancelled!
            externalCurrents(iteration) = externalCurrent/totalCurrent;
            totRECurrentDens(iteration) = currentConversionFactor*runawayCurrent; %In SI units
            totCurrentDens(iteration) = currentConversionFactor*totalCurrent; %In SI units
        end
    else
        totalCurrent = (xy2.*yWeights) * f1;

        if totalCurrent == 0
            currents(iteration) = 0;
            totRECurrentDens(iteration) = 0;
            totCurrentDens(iteration) = 0;
        else
            switch ss.runawayRegionMode
                case 0
                    runawayCurrent = (xy2.*yWeights.*nr_mask) * f1;
                case 1
                    runawayCurrent = CalculateRunawayMomentXiDependent(ss,fMinus1,xiIntMat,currentWeightsMat);
                otherwise
                    error('Invalid mode');
            end
            fastCurrent = (xy2.*yWeights.*tail_mask) * f1;
            currents(iteration) = runawayCurrent/totalCurrent; % Constants and normalization have been cancelled!
            fastCurrents(iteration) = fastCurrent/totalCurrent; % Constants and normalization have been cancelled!
            totRECurrentDens(iteration) = currentConversionFactor*runawayCurrent; %In SI units
            totCurrentDens(iteration) = currentConversionFactor*totalCurrent; %In SI units
        end
    end

    %Calculate the energy moment
    totalEnergy(iteration) = (yWeights.*energyMoment) * f0;
    fastEnergyFrac(iteration) = (yWeights.*energyMoment.*tail_mask) * f0 / totalEnergy(iteration);
    switch ss.runawayRegionMode
        case 0
            runawayEnergyFrac(iteration) = (yWeights.*energyMoment.*nr_mask) * f0 / totalEnergy(iteration);
        case 1
            en = CalculateRunawayMomentXiDependent(ss,fMinus1,xiIntMat,energyWeightsMat);
            runawayEnergyFrac(iteration) = en/totalEnergy(iteration);
        otherwise
            error('Invalid runaway region mode')
    end


    %%%%%% Prepare for matrix-building %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~sp.noTimeDep
        dt = sp.dts(iteration);
        veBar2 = sp.veBars(iteration)^2;
        veBar2NextStep = sp.veBars(min(iteration+1,end))^2;
        veBar3 = veBar2 * sp.veBars(iteration);
        nBar = sp.nBars(iteration);
        nBarNextStep = sp.nBars(min(iteration+1,end));
    end

    %(Re)build the matrix
    if  sp.parametersHaveChanged(iteration) || gridWasExtended
        %%% If parameters have changed, rebuild the matrix
        matrixRecomputations(iteration) = true;


        %Do the actual build
        if sp.noTimeDep
            %This distinction is needed when using automatic grid
            %extension to avoid extending the parameter vectors
            operator = BuildMatrix(ss,sp,1,x,y,dxdy,x2,y2,...
                ddy,d2dy2,ddx,d2dx2,deltaRef);
        else
            operator = BuildMatrix(ss,sp,iteration,x,y,dxdy,x2,y2,...
                ddy,d2dy2,ddx,d2dx2,deltaRef);
        end

        if includeSynchrotronLoss
            %The scaling is just a scalar
            if sp.noTimeDep
                operator = operator + synchrotronLossTerm;%+ (1/sp.nueeBars(1))*synchrotronLossTerm; A: %scaling inside object
            else
                operator = operator + synchrotronLossTerm;%(1/sp.nueeBars(iteration))*synchrotronLossTerm;
            end
        end

        switch ss.timeAdvanceMethod
            case 0
                % Backward Euler
                timeAdvanceMatrix = diagonalPlusBoundaryCondition - dt*operator;
                if ss.reduceMemoryConsumption
                    clear operator
                end
            case 1
                % BDF2
                timeAdvanceMatrix = diagonalPlusBoundaryCondition - (2/3)*dt*operator;
                if ss.reduceMemoryConsumption
                    clear operator
                end
            case 2
                % Trapezoid rule
                timeAdvanceMatrix = diagonalPlusBoundaryCondition - (dt/2)*operator;
            otherwise
                error('Invalid timeAdvanceMethod')
        end

        %Invert the matrix
        if ss.solveIteratively
            [factor_L,factor_U] = ilu(timeAdvanceMatrix,struct('type','ilutp','droptol',1e-1));
        else
            [factor_L, factor_U, factor_P, factor_Q] = lu(timeAdvanceMatrix);
            if ss.reduceMemoryConsumption
                clear timeAdvanceMatrix
            end
        end

    elseif sp.dtHasChanged(iteration) && iteration < nTimeSteps
        %If we use a logarithmic time step and the parameters have not
        %changed, we can save some effort by not rebuilding the operator,
        %just changing the time step and factorizing the matrix.

        switch ss.timeAdvanceMethod
            case 0
                % Backward Euler
                timeAdvanceMatrix = diagonalPlusBoundaryCondition - dt*operator;
            case 1
                % BDF2
                timeAdvanceMatrix = diagonalPlusBoundaryCondition - (2/3)*dt*operator;
            case 2
                % Trapezoid rule
                timeAdvanceMatrix = diagonalPlusBoundaryCondition - (dt/2)*operator;
            otherwise
                error('Invalid timeAdvanceMethod')
        end

        %Invert the matrix
        if ss.solveIteratively
            [factor_L,factor_U] = ilu(timeAdvanceMatrix,struct('type','ilutp','droptol',1e-1));
        else
            [factor_L, factor_U, factor_P, factor_Q] = lu(timeAdvanceMatrix);
        end
    end

    % Build the right-hand side
    switch ss.timeAdvanceMethod
        case 0
            rhs = fMinus1;
        case 1
            rhs = (4/3)*fMinus1 - (1/3)*fMinus2;
        case 2
            rhs = fMinus1 + (dt/2)*operator*fMinus1;
    end

    % At the moment, I treat the source term explicitly. This may
    % compromise the accuracy of the high-order time advance schemes,
    % so at some point I should think about whether there is a more accurate approach.
    if useBoltzmannBasedOperator

        %sourceVector = boltzmannMatrix*fMinus1; is done inside object
        rhs = rhs + sourcesPreFactor*dt*sourceVector;
    end
    switch this.eqSettings.sourceMode
        case 0
            % Nothing to do.
        case 1
            % Rosenbluth-Putvinskii-like source
            switch ss.initialDistribution
                case {0,2}
                    if ss.keepTrackOfExternalRunaways
                        nrEffective = nr+nExternal;
                    else
                        nrEffective = nr;
                    end
                case 1
                    % Provide some 'seed' runaways:
                    nrEffective = max([1, nr]);
            end
            sourceVector = sp.sourceTimeFactor(iteration)*nrEffective*RosenbluthSourceVector;
            rhs = rhs + sourcesPreFactor * dt * sourceVector;
        case 2
            % Rosenbluth-Putvinskii-like source, but using the number
            % of fast particles instead of the number of runaways
            if ss.keepTrackOfExternalRunaways
                nFastEffective = nFast+nExternal;
            else
                nFastEffective = nFast;
            end
            sourceVector = sp.sourceTimeFactor(iteration)*nFastEffective*RosenbluthSourceVector;
            rhs = rhs + sourcesPreFactor * dt * sourceVector;
        case 3
            %Chiu-Harvey-like source, y-int implementation
            intOnGrid = zeros(ss.Ny,1);
            idYForSource = idYFirstAboveCrit:idYMaxForCHSource;

            for L = 0:ss.Nxi-1
                indices = L*ss.Ny+(1:ss.Ny);

                %We need to perform one y integral per grid
                %point where the source is non-vanishing.
                %Restrict with indices to avoid unnecessary computation
                intOnGrid(idYForSource) = ...
                    cHSourceMatrix(idYForSource,idYMinCHInt:end,L+1) ...
                    * (cHIntVector(idYMinCHInt:end)...
                    .* f0(idYMinCHInt:end)')';
                sourceVector(indices) = (2*L+1)*nBar*cHSourcePreFactor.*intOnGrid;
                %cHSourcePreFactor is a vector on the grid, so we need
                %to include it for each mode, rather than as a
                %constant below.
            end
            rhs = rhs + sourcesPreFactor*dt*sourceVector;
        case {4,5}
            %
        otherwise
            error('Invalid sourceMode')
    end


    %Add a particle source to account for changes in density
    if (~sp.noTimeDep && sp.nHasChanged(iteration))
        %Add a particle source (delta function in time) to account for
        %the change in density
        changeFactor = sourcesPreFactor*(nBar/nBarNextStep-1);
        rhs(1:ss.Ny) = rhs(1:ss.Ny) + changeFactor * nBar/veBar3 ...
            * exp(-y2'/veBar2) .* (y2'/veBar2 - 5/2);
    end


    % Now step forward in time.
    if ss.solveIteratively
        [soln,convFlag] = gmres(timeAdvanceMatrix,...
            rhs,7,1e-6,20,factor_L,factor_U,fMinus1);

        %Alternative solvers:
        %                     [soln,convFlag] = bicgstab(timeAdvanceMatrix,rhs,1e-6,20,factor_L,factor_U,fMinus1);
        %                     gmres(A,b,restart,tol,maxit,M)

        if convFlag
            warning('The iterative solver failed to converge');
        end
    else
        % The next line is equivalent to 'soln = matrix \ rhs', but much faster:
        soln = factor_Q * (factor_U \ (factor_L \ (factor_P * rhs)));

    end


    %If we want a calculation of just the primary runaways as well,
    %repeat the time stepping with the source turned off
    if ss.calculateForOnlyPrimaries
        % Compute "number of primary runaways":
        f0Prim = fMinus1Prim(1:ss.Ny);
        nInstPrim = (yWeights.*y2) * f0Prim;
        switch ss.runawayRegionMode
            case 0
                nrPrim = (yWeights.*y2.*nr_mask) * f0Prim;
            case 1
                nrPrim = CalculateRunawayMomentXiDependent(ss,fMinus1Prim,xiIntMat,densityWeightsMat);
            otherwise
                error('Invalid mode');
        end
        nrPrim = nrPrim/nInstPrim;
        % Note: nrPrim is really the density of primary runaways
        % normalized by the instantaneous total electron density
        % (n_r/n). Some prefactors have been cancelled, so nInst can
        % not be used directly as a measure of the density
        nrsPrim(iteration) = nrPrim;

        switch ss.timeAdvanceMethod
            case 0
                rhsPrim = fMinus1Prim;
            case 1
                rhsPrim = (4/3)*fMinus1Prim - (1/3)*fMinus2Prim;
            case 2
                rhsPrim = fMinus1Prim + (dt/2)*operator*fMinus1Prim;
        end

        if ss.yMaxBoundaryCondition ~= 3
            rhsPrim((1:ss.Nxi)*ss.Ny) = 0;
        end

        if ss.solveIteratively
            [solnPrim,convFlag] = gmres(timeAdvanceMatrix,...
                rhsPrim,7,1e-6,20,factor_L,factor_U,fMinus1Prim);
            if convFlag
                warning('The iterative solver failed to converge');
            end
        else
            % The next line is equivalent to 'soln = matrix \ rhs', but much faster:
            solnPrim = factor_Q * (factor_U \ (factor_L \ (factor_P * rhsPrim)));
        end

        fMinus2Prim = fMinus1Prim;
        fMinus1Prim = solnPrim;
    end

    fMinus2 = fMinus1;
    fMinus1 = soln;
end
%The final time step is not saved, but we need to run through the loop to
%calculate all the moments etc of the next to final dist. (This can be
%avoided with a proper restructuring of the code.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of time-advance loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process and save the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For now, I'll just use this finite different to compute dn/dt:
dnrdts = (nrs-circshift(nrs,[0,1]))./(sp.dts); %The use of dt here implies that the units are per reference collision time!
dnrdts(1) = 0; %Why??? How does this work with starting from previous dist?

if ss.calculateForOnlyPrimaries
    dnrdtsPrim = (nrsPrim-circshift(nrsPrim,[0,1]))./sp.dts;
    dnrdtsPrim(1) = 0;
end

if ~ss.silentMode
    fprintf('The matrix was recomputed %d time(s) and the source %d time(s).\n',sum(matrixRecomputations),sum(sourceRecomputations))
    fprintf('Total elapsed time: %g\n',toc(startTime))
    fprintf('***************************************************************\n')
end

%%%%%%%% removed saving implementation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% End of main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for building the avalanche source: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [idTail,tail_mask] = DetermineFastParticleRegion(ss,x,y,fMinus1,...
            nBar,veBar3,deltaRef,nr_mask)
        switch ss.fastParticleDefinition
            case 0 %Do not calculate fast particles
                idTail = [];
            case 1 %Based on the beginning of the tail
                [fAtXi1,~] = SumLegModesAtXi1(fMinus1,ss.Ny,ss.Nxi);
                maxw = nBar/veBar3*exp(-y.^2/veBar2);
                quota = fAtXi1'./maxw;
                idTail = find(quota>ss.tailThreshold,1);
            case 2 %Based on relative speed
                idTail = find(x > ss.relativeSpeedThreshold,1);
            case 3 %Based on absolute speed
                idTail = find(deltaRef*x > ss.absoluteSpeedThreshold,1);
            case 4 %Combination
                idRE = find(nr_mask,1);
                idRel = find(x > ss.relativeSpeedThreshold,1);
                idAbs = find(deltaRef*x > ss.absoluteSpeedThreshold,1);
                idTail = min([idRE,idRel,idAbs]);
            otherwise
                error('Invalid fast particle definition selected.');
        end

        if idTail
            tail_mask = y > y(idTail);
        else
            tail_mask = zeros(size(y));
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions related to the grid:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function [lastUseful,usefulFrac] = FindLastUsedGridPoint(f,Ny,Nxi,usefulThreshold)
        [FAt1,~] = SumLegModesAtXi1(f,Ny,Nxi);
        lastUseful = find(FAt1 < usefulThreshold,1); %Look at the whole F (all modes) at xi=1
        if isempty(lastUseful)
            usefulFrac = 1;
        else
            usefulFrac = lastUseful/Ny;
        end
        %Check also where the tail "begins"

    end

%%%%%%%%%% used for nExternals
    function flowRate = GetFlowThroughGridBoundary(F0,F1,F2,y,y2,x,x2,ddy,...
            deltaRef,nBar,erfs,yBoundary,...
            includeSynchrotronLoss,sp,it)

        EHat = sp.EHats(it);
        BHat = sp.BHatRef;
        nueeBar = sp.nueeBars(it);
        veBar = sp.veBars(it);

        dF0dy = (ddy*F0)';
        yb = find(y >= yBoundary, 1);
        veBar2 = veBar*veBar;
        expx2 = exp(-x2/veBar2);
        %         erfs = erf(x/veBar);
        psi = (erfs-2/(veBar*sqrt(pi))*x.*expx2) * veBar2 ./ (2*x2);

        %Calculate growth rate
        eFieldTerm = EHat/3*y2.*F1';
        collisionalTerm = 3*sqrt(pi)/4 * nueeBar*veBar*veBar2 * y2.*psi.*(2/veBar2*F0' + dF0dy./x);
        if includeSynchrotronLoss
            FTilde = 2/3*(F0'-F2'/5);
            synchrotronTerm = BHat^2/nueeBar * y.*y2.*sqrt(1+deltaRef^2*y2).*FTilde;

            dndrsTot = -4/(sqrt(pi)*nBar) * (eFieldTerm+collisionalTerm+synchrotronTerm); %This is (dn_r/dt) / (nu_ee n)
        else
            dndrsTot = -4/(sqrt(pi)*nBar) * (eFieldTerm+collisionalTerm); %This is (dn_r/dt) / (nu_ee n)
        end

        flowRate = dndrsTot(yb);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for the xi-dependent runaway region:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [m,xiGrid] = GenerateCumIntMat(Nxi,nPoints)
        %Cumulatively integrate legendre modes over a xi grid

        xiGrid = linspace(0,1,nPoints);
        dXi = xiGrid(2)-xiGrid(1);
        legPols = LegendrePolynomials(Nxi-1,xiGrid);
        m = dXi*fliplr(cumtrapz(fliplr(legPols),2));
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for screening:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = MapToXiIntMat(cumintMat,xiGrid,EOverEc,y2,deltaRef)
    %Assumes that xiGrid is uniformly spaced, starting at 0!

    xics = 1/EOverEc * (1./(deltaRef^2*y2) + 1);
    xics(xics>1) = 1; %Remove points close to the bulk that give unphysical results. They are not in the runaway region, anyway, and will not contribute in the end
    dxi = xiGrid(2)-xiGrid(1);
    xiIds = round(xics/dxi)+1; %Find the id of the nearest point (assumes uniformly spaced xiGrid starting at 0!)
    m = cumintMat(:,xiIds);


    %         figure(888);
    %         y = sqrt(y2);
    %         yPara = xics.*y;
    %         yPerp = y.*sqrt(1-xics.*xics);
    %         plot(yPara,yPerp);
    %         xlabel('y_{||}')
    %         ylabel('y_\perp')
    %         legend('y_c')
    %         drawnow
    %         pause(0.1)

    %Each column corresponds to a specific y value. Set those for y<y_c
    %to 0 so that they don't contribute to the integral
    yc = 1/(deltaRef*sqrt(EOverEc-1));
    id = find(y2>yc^2,1);
    m(:,1:id) = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end
