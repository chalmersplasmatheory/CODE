%*******************************************************************************************************************from raw run function

%Print info about the run
if ~ss.silentMode
    fprintf('***************************************************************\n')
    fprintf('Time-dependent CODE, using ')
    switch ss.timeAdvanceMethod
        case 0
            % Backward Euler
            fprintf('backward Euler')
        case 1
            % BDF2
            fprintf('BDF2')
        case 2
            % Trapezoid rule
            fprintf('trapezoid rule')
        otherwise
            error('Invalid timeAdvanceMethod')
    end
    fprintf(' for time advance.\n')

    fprintf('Using the ');
    switch ss.collisionOperator
        case 0
            fprintf('relativistic');
        case 1
            fprintf('non-relativistic, momentum conserving');
        case 2
            fprintf('combined approximate relativistic, momentum conserving');
        case 3
            fprintf('approximate ad-hoc relativistic, momentum conserving');
        otherwise
            error('Invalid collision opeartor.');
    end
    fprintf(' collision operator.\n');

    switch ss.sourceMode
        case 0
            fprintf('No avalanche source is used.\n');
        case {1,2}
            fprintf('A Rosenbluth-Putvinski-like avalanche source is used.\n');

            if ss.sourceMode == 2
                fprintf('The avalanche source is determined by the population of fast particles,\ndefined ');
            elseif ss.fastParticleDefinition
                fprintf('Fast particles are defined ');
            end
            switch ss.fastParticleDefinition
                case 0
                    %Don't print anything
                case 1
                    fprintf('by the start of the ''tail''.\n');
                case 2
                    fprintf('by the relative particle speed (v/v_e > %.1f)\n',...
                        ss.relativeSpeedThreshold);
                case 3
                    fprintf('by the absolute particle speed (v/c > %.1f)\n',...
                        ss.absoluteSpeedThreshold);
                case 4
                    fprintf('using the most generous of v>v_crit, v/v_e > %.1f and v/c > %.1f.\n',...
                        ss.relativeSpeedThreshold,...
                        ss.absoluteSpeedThreshold);
                otherwise
                    error('Invalid fast particle definition.\n')
            end
        case {3,4,5}
            fprintf('A cross-section dependent avalanche source is used. ');
            switch ss.sourceMode
                case 3
                    fprintf('Only F_0 of the RE distribution is used (Chiu-Harvey).\n');
                case 4
                    fprintf('The complete RE distribution is used (Embreus), but no sinks are included.\n');
                case 5
                    fprintf('The complete RE distribution is used (Embreus), and the sinks are included.\n');
            end
        otherwise
            error('Invalid avalanche source mode.');
    end

    if exist('gridParam','var')
        fprintf('Using a non-uniform grid (mode %d) with grid parameter %.2g\n',ss.yGridMode,ss.gridParameter);
    end

end

deltaRef = ComputeDeltaFromT(ss.TRef);


%%%%% 'Live' plots %%%%%
%Make a plot for diagnosing the time dependence of the parameters
if ss.showVarPlot
    if sp.allScalar %Only plot the time-dependent parameters if they are not constant!

        ss.showVarPlot = false;
        if ishandle(2 + ss.figureOffset)
            close(2 + ss.figureOffset);
            %Close an already open figure, to avoid confusion in showing
            %outdated data
        end

    else


        figure(2 + ss.figureOffset);
        clf;

        nRowsVar = 5;
        nColsVar = 2;
        finalPos = [9,10];
        hProgressLines = zeros(1,9);
        progLineColor = '-k';
        t0 = sp.timestepsForPlot(1);

        if nargin > 0 && os.timeStepMode<2
            %Read in parameter timevectors again, since they have been altered
            timestepsT = os.tT;
            timestepsN = os.tn;
            timestepsZ = os.tZ;
            timestepsE = os.tE;
        end

        [inTemps, timestepsT] = ExpandForPlot(os.T,timestepsT,sp.timestepsForPlot);
        [inDens, timestepsN] = ExpandForPlot(os.n,timestepsN,sp.timestepsForPlot);
        [inEFields, timestepsE] = ExpandForPlot(os.E,timestepsE,sp.timestepsForPlot);
        [inZs, timestepsZ] = ExpandForPlot(os.Z,timestepsZ,sp.timestepsForPlot);

        %%%Temp, n, time step related plots
        %Temp
        subplot(nRowsVar,nColsVar,1)
        %Make a linear or semilog plot, depending on the range of the input
        tChange = max(sp.temperatures)/min(sp.temperatures);
        if tChange < 20
            hold on;
            plot(timestepsT,inTemps,'-k','linewidth',2);
            plot(sp.timestepsForPlot,sp.temperatures,'+g');
            plot(sp.timestepsForPlot([1,end]),[ss.TRef,ss.TRef],'--b');
        else
            semilogy(timestepsT,inTemps,'-k','linewidth',2);
            hold on;
            semilogy(sp.timestepsForPlot,sp.temperatures,'+g');
            semilogy(sp.timestepsForPlot([1,end]),[ss.TRef,ss.TRef],'--b');
            range = round(log10(0.5*min(sp.temperatures))):round(log10(5*max(sp.temperatures)));
            set(gca,'ytick',10.^range);
        end
        set(gca,'ylim',[0.9,1.08].*get(gca,'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        legend('Input','Interpolated','Reference (grid)','location','NorthEast');
        ylabel('Temperature [eV]')
        hProgressLines(1) = plot([t0,t0],get(gca,'YLim'),progLineColor);
        %Density
        subplot(nRowsVar,nColsVar,3)
        hold on;
        plot(timestepsN,inDens,'-k','linewidth',2);
        plot(sp.timestepsForPlot,sp.densities,'+r');
        plot(sp.timestepsForPlot([1,end]),[ss.nRef,ss.nRef],'--b');
        set(gca,'ylim',[0.9,1.08].*get(gca,'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        legend('Input','Interpolated','Reference (grid)','location','Best');
        ylabel('Density [m^{-3}]')
        hProgressLines(2) = plot([t0,t0],get(gca,'YLim'),progLineColor);
        %Thermal speed
        subplot(nRowsVar,nColsVar,5)
        hold on;
        plot(sp.timestepsForPlot,sp.veBars*deltaRef,'-','Color',[1,0.8,0]);
        plot(sp.timestepsForPlot([1,end]),[deltaRef,deltaRef],'--b');
        lims = get(gca,'ylim');%
        set(gca,'ylim',[0.9,1.08].*get(gca,'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        legend('Instantaneous','Reference','location','best');
        ylabel('Thermal speed (v_e/c)')
        hProgressLines(3) = plot([t0,t0],lims,progLineColor);
        %Time step
        subplot(nRowsVar,nColsVar,7)
        hold on;
        [ax,h1,h2] = plotyy(sp.timestepsForPlot,sp.dts,sp.timestepsForPlot,sp.dts.*sp.nueeBars);
        title('Time step');
        set(get(ax(1),'Ylabel'),'String','in reference coll. times')
        set(get(ax(2),'Ylabel'),'String','in instantaneous coll. times')
        set(h1,'Linestyle','--','linewidth',1);
        set(h2,'Linestyle','-','linewidth',2);
        set(ax(1),'ylim',[0.9,1.08].*get(ax(1),'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        set(ax(2),'ylim',[0.9,1.08].*get(ax(2),'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        xlabel(sp.timeLabelText);

        hProgressLines(4) = plot([t0,t0],get(gca,'YLim'),progLineColor);



        %%% E, Z, lnLambda
        %E
        subplot(nRowsVar,nColsVar,2)
        hold on;
        plot(timestepsE,inEFields,'-k','linewidth',2);
        plot(sp.timestepsForPlot,EFields,'+m');
        legend('Input','Interpolated','location','NorthEast');
        ylabel('Electric field [V/m]')
        set(gca,'ylim',[0.9,1.08].*get(gca,'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        hProgressLines(5) = plot([t0,t0],get(gca,'YLim'),progLineColor);
        %E/E_c, EHat
        subplot(nRowsVar,nColsVar,4)
        hold on;
        [ax,h1,h2] = plotyy(sp.timestepsForPlot,sp.EOverEc,sp.timestepsForPlot,sp.EHats);
        %             legend([h1,h2],{'E/E_c','EHat'},'location','NorthEast');
        set(get(ax(1),'Ylabel'),'String','E/E_c')
        set(get(ax(2),'Ylabel'),'String','EHat \propto E/E_D')
        set(h1,'Linestyle','-','linewidth',2);
        set(h2,'Linestyle','--','linewidth',2);
        set(ax(1),'ylim',[0.9,1.08].*get(ax(1),'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        set(ax(2),'ylim',[0.9,1.08].*get(ax(2),'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        hProgressLines(6) = plot([t0,t0],get(gca,'YLim'),progLineColor);
        %Z
        subplot(nRowsVar,nColsVar,6)
        hold on;
        plot(timestepsZ,inZs,'-k','linewidth',2);
        plot(sp.timestepsForPlot,sp.Zs,'+c');
        legend('Input','Interpolated','location','NorthEast');
        ylabel('Z')
        set(gca,'ylim',[0.9,1.08].*get(gca,'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        hProgressLines(7) = plot([t0,t0],get(gca,'YLim'),progLineColor);
        %lnLambda
        subplot(nRowsVar,nColsVar,8)
        %             [~,~,~,~,lnLambdas,~,~] = DerivedParameters(1,sp.temperatures,sp.densities,[]);
        [~,lnLambdaRef,~,~] = getDerivedParameters(ss.TRef,ss.nRef,[]);
        hold on;
        plot(sp.timestepsForPlot,sp.lnLambdas,'+m');
        plot(sp.timestepsForPlot([1,end]),[lnLambdaRef,lnLambdaRef],'--b');
        set(gca,'ylim',[0.9,1.08].*get(gca,'ylim'),'xlim',sp.timestepsForPlot([1,end])); %Make some vertical room
        legend('Interpolated','Reference (grid)','location','NorthEast');
        xlabel(sp.timeLabelText);
        ylabel('ln\Lambda')
        hProgressLines(8) = plot([t0,t0],get(gca,'YLim'),progLineColor);



        stringForTop='Time variation of plasma parameters and time steps';
        annotation('textbox',[0 0.95 1 .08],'HorizontalAlignment',...
            'center','Interpreter','none','VerticalAlignment','bottom',...
            'FontSize',11,'LineStyle','none','String',stringForTop);




        % Progress bar / matrix recomputations
        hProgressAx = subplot(nRowsVar,nColsVar,finalPos);
        hold on;
        progressAxLim = [0,1];
        set(hProgressAx,'XLim',[t0,sp.timestepsForPlot(end)]);
        set(hProgressAx,'YLim',progressAxLim);
        progAxPos = get(hProgressAx,'Position');
        set(hProgressAx,'Position',[progAxPos(1:3),progAxPos(4)/2]);
        for i = 1:length(sp.timestepsForPlot);
            plot([sp.timestepsForPlot(i),sp.timestepsForPlot(i)],progressAxLim,'--r');
        end
        hProgressLines(end) = plot([t0,t0],progressAxLim,'-b');
        xlabel(sp.timeLabelText)
        title('Time steps (red), together with matrix (black) and source (yellow) recomputations and grid extensions (green)')

    end
end





% *******************************************************************************************************************************from ProcessInputs



%             timesteps = 0:this.plasma.dt:this.plasma.tMax;% This can't be right, right?
% Create a timestepsForPlot vector
switch this.plasma.timeUnit
    case 'normalized'
        timestepsForPlot = timesteps;
        timeLabelText = 'time (1/\nu_{ee,ref})';
    case 's'
        timestepsForPlot = timesteps/nueeRef;
        timeLabelText = 'time (s)';
    case 'ms'
        timestepsForPlot = 1e3*timesteps/nueeRef;
        timeLabelText = 'time (ms)';
end



if this.plasma.showRunawayRatePlot
    %Some parameter vectors need to have one entry per time
    %step if we want to show the runaway rate plot
    requiredLength = length(timesteps);
else
    requiredLength = 2;
end


% **********************************************************************************************************************from directly after ProcessInputs

%Initiate a plot to show the time evolution of the distribution.
if ss.showDistEvolutionPlot
    hDistEvo = figure(3 + ss.figureOffset);
    clf(hDistEvo);
    hDistEvoAxes = axes;
    xlim([-round(0.3*ss.yMaxForPlot),ss.yMaxForPlot])
end

%Initialize a plot to show the avalanche source term
if ss.showAvalancheSourcePlot
    hAvaSourceFig = figure(13 + ss.figureOffset);
    clf(hAvaSourceFig);
end

%Determine when to update the "live" plots
isUpdateTimestep = zeros(size(sp.timesteps));
if ss.showVarPlot || ss.showDistEvolutionPlot || ss.showAvalancheSourcePlot
    isUpdateTimestep(2:ss.stepSkip:end) = 1;
end

%********************************************************************************************************************from row 723 and onwards in actual CODE
%%% If an initial distribution is provided, rescale it if necessary.
%Otherwise initialize a grid from scratch
usefulFrac = 0;
gridIsInitialized = 0;
initialStepDensityChange = 0;
initialStepTempChange = 0;
if this.plasma.initialDistribution == 2

    sameNy = (inDist.Ny == this.plasma.Ny);
    sameYMax = (inDist.yMax == this.plasma.yMax);
    sameNxi = (inDist.Nxi == this.plasma.Nxi);
    sameTRef = (inDist.TRef == this.plasma.TRef);
    sameNRef = (inDist.nRef == this.plasma.nRef);
    sameYGridMode = (inDist.yGridMode == this.plasma.yGridMode);
    % if gridMode == 5, we need all grid parameters to be the same. Note
    % that here gridStepPosition is normalized to yc, which depends on Ec
    sameParam = (inDist.settings.gridParameter == this.plasma.gridParameter);
    sameWidth = (inDist.settings.gridStepWidth == this.plasma.gridStepWidth);
    samePos = (inDist.settings.gridStepPosition == this.plasma.gridStepPosition) && ...
        (inDist.EOverEc(1) == sp.EOverEc(1));
    sameGridParameter = (this.plasma.yGridMode~= 5) || (sameParam && sameWidth && samePos);

    %Rescale the grid if it doesn't match the current one
    if sameNy && sameYMax && sameNxi && sameTRef && sameYGridMode && sameGridParameter
        %Does match

        providedInitialDist = inDist.f(:,step);
        gridNormGridPoint = inDist.gridNormGridPoint;
        gridNormVal = inDist.gridNormVal;

        %Check if the tail has reached the end of the grid
        if this.plasma.useAutomaticGridExtension
            [doExtend,usefulFrac] = CheckForGridExtension(providedInitialDist,this.plasma.Ny,this.plasma.Nxi,this.plasma.usefulThreshold);
            if doExtend
                [providedInitialDist,y,ddy,d2dy2,yWeights,this.plasma,matrixSize,~] = ...
                    DoExtendGrid(providedInitialDist,this.plasma,gridNormGridPoint,...
                    gridNormVal,sp.maxYForGridGrowth,0);
                gridIsInitialized = 1;
                this.plasma.yMaxForPlot = this.plasma.yMax;
                if this.plasma.updateSettingsObject
                    settings.Nxi = this.plasma.Nxi;
                    settings.Ny = this.plasma.Ny;
                    settings.yMax = this.plasma.yMax;
                    settings.yMaxForPlot = this.plasma.yMaxForPlot;
                end
            end
        end
        if ~sameNRef
            %Only nRef has changed, and there is no need to remake the
            %grid. Just rescaling F is enough
            providedInitialDist = inDist.nRef/this.plasma.nRef*providedInitialDist;
        end




        if ~this.plasma.silentMode
            fprintf('An initial distribution was provided.\n');
        end

    else
        %We need to redo everything
        fprintf('Reshaping the grid.');

        dist = inDist.f(:,step);
        inNxi = inDist.Nxi;
        inNy = inDist.Ny;
        inyMax = inDist.yMax;
        iny = inDist.y';
        inTRef = inDist.TRef;
        inNRef = inDist.nRef;

        if ~sameTRef
            %When decreasing TRef, make sure we keep all the
            %useful information in f. This might mean that we have to
            %change yMax/Ny/Nxi.
            fprintf('\b due to a change in temperature. ');
            yRatio = this.plasma.Ny/this.plasma.yMax;
            oldNy = this.plasma.Ny;
            [lastUseful,~] = FindLastUsedGridPoint(dist,inNy,inNxi,this.plasma.usefulThreshold);
            if lastUseful
                if inNy == lastUseful
                    yLast = iny(end);
                else
                    yLast = iny(lastUseful+1);
                end
                if this.plasma.TRef<inTRef
                    %@@@This is only good if the tail has stopped
                    %growing, but this is not always the case. Add a
                    %condition on EHAT to make sure?
                    newYMax = max([this.plasma.yMax,round(1.3*sqrt(inTRef/this.plasma.TRef)*yLast)]); %Make sure the new grid includes the info in f_0 (and an extra 30%), but not more
                else
                    newYMax = max([this.plasma.yMax,round(2*sqrt(inTRef/this.plasma.TRef)*yLast)]); %We need some extra grid space when increasing T @@@(?)
                end
            else
                newYMax = round(sqrt(inTRef/this.plasma.TRef)*this.plasma.yMax); %Otherwise use all the points
            end
            %Make sure that we don't make the grid larger than allowed by
            %externalBoundaryMeV
            newYMax = min(newYMax,sp.maxYForGridGrowth);


            switch this.plasma.yGridMode
                case 4
                    %Make use of the specific non-uniformity of the
                    %grid (y=s^2+gridParameter*s) to keep the
                    %resolution with the minimum amount of grid points.
                    NyNorm = inDist.gridNormGridPoint; %Here we need to use the grid normalization point and value to get the proper result for dy
                    yMaxNorm = inDist.gridNormVal;
                    c1 = yMaxNorm/newYMax * (1/NyNorm^2 + this.plasma.gridParameter/NyNorm);
                    this.plasma.Ny = ceil(1/c1 * (0.5*this.plasma.gridParameter + sqrt(0.25*this.plasma.gridParameter^2+c1))); %Scale the number of grid points
                    this.plasma.yMax = newYMax;
                    this.plasma.Nxi = ceil(this.plasma.Nxi+this.plasma.extensionXiSlope*(this.plasma.Ny-oldNy));
                otherwise
                    %Just scale the number of points linearly
                    this.plasma.Ny = round(newYMax*yRatio); %Scale the number of grid points
                    this.plasma.Nxi = round(this.plasma.Nxi+this.plasma.extensionXiSlope*(this.plasma.Ny-oldNy));
                    this.plasma.yMax = newYMax;
            end

            if this.plasma.updateSettingsObject
                settings.Nxi = this.plasma.Nxi;
                settings.Ny = this.plasma.Ny;
                settings.yMax = this.plasma.yMax;
                settings.yMaxForPlot = this.plasma.yMax;
            end

            fprintf('New size: Ny=%d, yMax=%d, Nxi=%d\n',this.plasma.Ny,this.plasma.yMax,this.plasma.Nxi);
        end

        %If TRef changes, F must be renormalized
        fRescaleFactor = inNRef/this.plasma.nRef * (this.plasma.TRef/inTRef)^1.5;

        [y,ddy,d2dy2,yWeights] = InitializeGrid(this.plasma.yGridMode,this.plasma.Ny,this.plasma.yMax,...
            this.plasma.gridParameter,this.plasma.gridStepWidth,this.plasma.gridStepPosition);
        gridNormGridPoint = ss.Ny;
        gridNormVal = ss.yMax;
        gridIsInitialized = 1;
        matrixSize = ss.Nxi*ss.Ny;
        providedInitialDist = zeros(matrixSize,1);
        if ss.showGridInterpolationPlot
            figure(8 + ss.figureOffset);
            clf;
            %Top text
            annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
                'Interpreter','none','VerticalAlignment','bottom',...
                'FontSize',11,'LineStyle','none','String',...
                'Initial distribution represented on native and rescaled grids');
        end
        %For each Legendre mode, interpolate the provided F_l to
        %the current grid. If the current run has more Legendre modes
        %than the initial one, leave them empty.
        for l = 0:min(inNxi,ss.Nxi)-1
            distrL = dist(l*inNy+(1:inNy));
            if ss.showGridInterpolationPlot && l<=3
                subH = subplot(2,2,l+1);
                %When plotting we rescale the old grid and F, so that we can
                %look at overlapping distributions to judge the
                %quality of the interpolation.
                %(This applies when inTRef ~= TRef || inNREf~=ss.B)
                semilogy(subH,sqrt(inTRef/ss.TRef)*iny,fRescaleFactor*distrL,'b');
                hold on;
                pH = semilogy(subH,sqrt(inTRef/ss.TRef)*iny, -fRescaleFactor*distrL, ':b') ;
                set(get(get(pH,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off'); % Exclude line from legend
            end
            if ss.yMax > sqrt(inTRef/ss.TRef)*inyMax
                %If the current grid extends further in y than the
                %provided initial grid, add a point corresponding to
                %the maximum of the current grid to the initial grid,
                %so that we can interpolate without problems
                distrL = fRescaleFactor * interp1([sqrt(inTRef/ss.TRef)*iny;ss.yMax],...
                    [distrL;distrL(end)],y,'pchip');
            else
                distrL = fRescaleFactor * interp1(sqrt(inTRef/ss.TRef)*iny,...
                    distrL,y,'pchip');
            end
            %Build the initial distribution on the current grid
            providedInitialDist(l*ss.Ny+(1:ss.Ny)) = distrL;
            if ss.showGridInterpolationPlot && l<=3
                semilogy(subH,y,distrL,'r');
                pH = semilogy(subH,y, -distrL, ':r');
                set(get(get(pH,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off'); % Exclude line from legend
                hold off;
                set(subH,'XLim',[y(1),max(iny(end),y(end))]);
                xlabel('y [in units of the new grid]');
                ylabel('F_l');
                title(['Legendre mode l = ',num2str(l)]);
                legend(subH,'On native grid','On grid used in this run',...
                    'Location','SouthWest');
            end
        end
        if ss.showGridInterpolationPlot
            drawnow
        end
        if ~ss.silentMode
            fprintf('The provided initial distribution was remapped to the current grid.\n')
        end
    end

    if ss.showDistEvolutionPlot
        %Change the plot limits in case the grid has changed
        xlim(hDistEvoAxes,[-round(0.3*ss.yMaxForPlot),ss.yMaxForPlot])
    end

    %Take care of the case where we do a restart and the density or
    %temperature has changed compared to the initial distribution.
    if inDist.settings.n(end) ~= ss.n(1)
        initialStepDensityChange = 1;
    end
    tempOfInitialDist = inDist.veBars(min(step,end))*inDist.delta/deltaRef; %Check also with different TRef to be sure
    if tempOfInitialDist ~= sp.veBars(1)
        initialStepTempChange = 1;
    end
else
    gridNormGridPoint = ss.Ny;
    gridNormVal = ss.yMax;
end

%******************************************************************* taken from somewhere a bit below ending of last

if ss.showMeVPlot
    figure(11+ss.figureOffset);
    clf;
    mev = 0.511 * (sqrt(1+deltaRef^2*y.^2)-1);
    plot(y,mev,'-ok');
    grid on
    xlabel('y');
    ylabel('E_k [MeV]');
end





%********************************************************************* taken from main loop end, where all plots are, basically between the two comments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


colors = [    0,   0, 1;
    1,   0, 0;
    0, 0.8, 0;
    0.8, 0.8, 0;
    0.5, 0.5, 0.5;
    1,   0, 1];

if ss.showGrids % Useful when working with non-uniform grids

    figure(7+ss.figureOffset)
    clf
    runsToPlot = 1;
    legendText={};

    for i = runsToPlot
        plot(y, -i*ones(size(y)),'.','Color',colors(i,:))
        hold on
    end
    xlabel('y')
    if numel(legendText)>0
        ylabel('-run number')
        legend(legendText)
    end
    ylim([-max(runsToPlot)-1,0])
end


if ss.showSolution %The distribution vector in the last time-step
    figure(5 + ss.figureOffset)
    clf
    semilogy(abs(soln))
    title('solution vector')
    xlabel('row in linear system')
end


if ss.showMainPlot

    figure(1 + ss.figureOffset)
    clf
    numRows=2;
    numCols=3;
    numRuns = 1;

    warning off MATLAB:Axes:NegativeDataInLogAxis

    shortLegendText={'Base case'};
    longLegendText = shortLegendText;
    longLegendText{end+1} =  'exp(-y^2)';

    for L=0:2
        subplot(numRows,numCols,L+1)
        for runNum=1:numRuns
            vectorToPlot = soln(L*ss.Ny + (1:ss.Ny));
            semilogy(y', vectorToPlot, 'Color',colors(mod(runNum-1,size(colors,1))+1,:))
            hold on
            plotHandle = semilogy(y', -vectorToPlot, ':', 'Color',colors(mod(runNum-1,size(colors,1))+1,:));
            set(get(get(plotHandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
        legendText = shortLegendText;
        if L==0
            semilogy(y, exp(-y.^2), '--k')
            legendText = longLegendText;
        end
        ylim([ss.FMinForPlot, 1])
        title(['L = ',num2str(L),' Legendre mode of distribution'])
        legend(legendText)
        xlabel('y = \gamma v / v_{th}')
        ylabel('F')

    end

    subplot(numRows,numCols,[4 5])
    for runNum=1:numRuns
        %Positive and negative here has been reversed since the code was
        %written. Negative X here thus means the direction of the runaway
        %tail.
        FPositiveX = zeros(ss.Ny,1);
        FNegativeX = zeros(ss.Ny,1);
        for L=0:(ss.Nxi-1)
            vectorToAdd = soln(L*ss.Ny + (1:ss.Ny));
            FPositiveX = FPositiveX + vectorToAdd;
            if mod(L,2)==0
                FNegativeX = FNegativeX + vectorToAdd;
            else
                FNegativeX = FNegativeX - vectorToAdd;
            end
        end
        xBig = -[fliplr(-y), y]';
        F = [flipud(FNegativeX); FPositiveX];
        semilogy(xBig, F, 'Color',colors(mod(runNum-1,size(colors,1))+1,:))
        hold on
        plotHandle = semilogy(xBig, -F, ':','Color',colors(mod(runNum-1,size(colors,1))+1,:));
        set(get(get(plotHandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
    end
    semilogy(xBig, exp(-xBig.^2), '--k')
    title('Total distribution for v_{perp} = 0')
    legend(longLegendText)
    xlabel('\gamma v_{||} / v_{th}')
    ylabel('F')
    %     xlim([-round(0.1*ss.yMaxForPlot),ss.yMaxForPlot])
    xlim([-10,ss.yMaxForPlot])
    ylim([ss.FMinForPlot, 1])


    subplot(numRows,numCols,6)

    %Compare with Spitzer resistivity (for the first run)
    if sp.noTimeDep
        resTemps = this.plasma.T(1);
        resNs = this.plasma.n(1);
        resEOEcs = sp.EOverEc(1);
        resZs = this.plasma.Z(1);
    else
        resTemps = this.plasma.T;
        resNs = this.plasma.n;
        resEOEcs = sp.EOverEc;
        resZs = this.plasma.Z;
    end
    [~,resLnLambdas,~,~] = getDerivedParameters(resTemps,resNs,[]);
    [oneOverEc,~, ~] = getNormalizedEFields(1,resTemps,resNs);
    Es = resEOEcs./oneOverEc; %We can't use the input E directly since it is defined at different times
    m0TimesE = 9.109380e-31 * 1.602176e-19; %SI - for electron
    eps0Sq = (8.854188e-12)^2; %SI
    etaNormWOZs = 4*sqrt(2*pi*m0TimesE)*resLnLambdas./(16*pi^2*eps0Sq*3*(resTemps).^(3/2));
    etaNormWZs = 4*sqrt(2*pi*m0TimesE*resZs).*resLnLambdas./(16*pi^2*eps0Sq*3*(resTemps).^(3/2));
    etaHatWZs = Es./(totCurrentDens.*etaNormWZs);
    etaHatWOZs = Es./(totCurrentDens.*etaNormWOZs);

    hold on
    plot(sp.timesteps(5:end), etaHatWOZs(5:end), 'Color','b')
    plot(sp.timesteps(5:end), etaHatWZs(5:end),'-.', 'Color',[0,0.5020,0])
    xLim = get(gca,'xlim');
    ylim([0,2.5]);
    plot(xLim,[0.99651105,0.99651105],'--','Color',[1,0.2,0.2]);
    plot(xLim,[0.50611832, 0.50611832],'--','Color',[0,0,0]);
    legend('Without Z^{1/2} factor','With Z^{1/2} factor','Value with model op.','Spitzer resistivity','Location','best')
    title('Normalized resistivity')
    xlabel('time (1/\nu_{ee,ref})')






    stringForTop='Results at the end of simulation(s).  Solid curves are positive, dashed curves are negative.';
    annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',11,'LineStyle','none','String',stringForTop);

end


if ss.showContourPlots %Contour plot of the distribution (final time-step)

    figure(6 + ss.figureOffset)
    clf
    numContours = 15;

    % Maximum speed which can fit in a square contour plot without having to
    % extrapolate, rounding sqrt(2) up slightly:
    yMax2 = ss.yMax/1.5;

    %     yMax2s = [yMax2, yMax2*0.5, yMax2*0.2];
    yMax2s = yMax2*0.2;

    for plotNum = 1:numel(yMax2s)

        % Number of grid points to use in making the contour plots:
        NPlot = 81;
        % It's best to choose NPlot odd to avoid the singularity at the origin.

        yMax2 = yMax2s(plotNum);
        %         yPar1D = linspace(-yMax2, yMax2, NPlot);
        %         yPerp1D = yPar1D;
        yPar1D = linspace(-2.5*yMax2, 0.5*yMax2, NPlot);
        yPerp1D =linspace(-yMax2, yMax2, NPlot);
        [yPar2D, yPerp2D] = meshgrid(yPar1D, yPerp1D);
        y2D = sqrt(yPar2D.*yPar2D + yPerp2D.*yPerp2D);
        xi2D = yPar2D ./ y2D;

        y2D = reshape(y2D,[NPlot*NPlot,1]);
        xi2D = reshape(xi2D,[NPlot*NPlot,1]);

        f = zeros(NPlot*NPlot,1);
        for L=0:(ss.Nxi-1)
            fSlice = soln(L*ss.Ny+(1:ss.Ny));

            % Use Matlab's built-in function to evaluate the Legendre
            % polynomials on the points we need for the contour plot:
            LegendreMatrix = legendre(L, xi2D);
            % Matlab's 'legendre' function also returns a bunch of other
            % uninteresting numbers besides the Legendre polynomials, so
            % keep only the values of the Legendre polynomials:
            Legendres = LegendreMatrix(1,:)';

            % Interpolate from the regular y grid to the points we need in
            % the Cartesian plane for the contour plot:
            f = f + Legendres .* interp1(y, fSlice, y2D);
        end

        f = reshape(f, [NPlot, NPlot]);

        toPlot = log10(abs(f));
        minContour = max([-15, min(min(toPlot))]);
        contours = linspace(minContour, 0, numContours);%linspace(minContour, max(max(toPlot)), numContours);

        %         subplot(2,2,plotNum)

        contourf(-yPar2D, yPerp2D, toPlot, contours)
        colorbar
        xlabel('y_{||}')
        ylabel('y_{perp}')
        %         xlabel('\gamma v_{||} / v_{th}')
        %         ylabel('\gamma v_{perp} / v_{th}')
        title('log_{10}|F|')


        hold on;
        %Isotropic
        yc = wMinForSource/deltaRef;
        th = linspace(0,2*pi,1000);
        circPar = cos(th);
        circPerp = sin(th);
        h1 = plot(yc*circPar,yc*circPerp,'--r','LineWidth',2);

        %Region in Smith et al. (valid non-realtivistically)
        %         xiSep = 2*xc2./x2 - 1;
        xi = circPar;
        EED = 2/(3*sqrt(pi))*sp.EHats(end);
        xc2 = 0.5/EED;
        xSep = sqrt(2*xc2./(xi+1));
        xSPar = xi.*xSep;
        xSPerp = xSep.*sqrt(1-xi.*xi);
        xSPerp(xSPar<-10) = [];
        xSPar(xSPar<-10) = [];
        xSPar (xSPerp>max(yPerp2D(:))) = [];
        xSPerp(xSPerp>max(yPerp2D(:))) = [];
        plot(xSPar,xSPerp,'-.c','LineWidth',2);
        h3 = plot(xSPar,-xSPerp,'-.c','LineWidth',2);

        %Xi-dependent
        xi = circPar(circPar>0);
        arg = xi*sp.EOverEc(end)-1;
        xi(arg<=0) = [];
        arg(arg<=0) = [];
        ycX = 1./(deltaRef*sqrt(arg));
        ycPar = xi.*ycX;
        ycPerp = ycX.*sqrt(1-xi.*xi);
        ycPar(ycPerp>max(yPerp2D(:))) = [];
        ycPerp(ycPerp>max(yPerp2D(:))) = [];
        plot(ycPar,ycPerp,'--y','LineWidth',2);
        h2 = plot(ycPar,-ycPerp,'--y','LineWidth',2);



        legend([h1,h2,h3],'y_c (isotropic)','y_c(\xi)','x_{sep} (Smith)')

        axis equal
    end
end


if ss.showRunawayRatePlot %Plot showing runaway generation over time.
    %New plot showing runaway generation over time. Should correlate with
    %E-field strength. Only plots the data for the first run, if a scan is
    %performed.

    %Estimate the generation from analytical formulas (per reference collision time )
    dreicerEstimate = EstimateDreicerGeneration(sp.EOverED,this.plasma.Z,sp.nueeBars);
    avalancheEstimate = EstimateAvalancheGeneration(sp.EOverEc,nrs,...
        sp.lnLambdas,sp.deltas,this.plasma.Z,sp.nueeBars);
    %     dreicerEstimate = EstimateDreicerGeneration(op.EOverED(1:end-1),op.Zs(1:end-1));
    %     avalancheEstimate = EstimateAvalancheGeneration(op.EOverEc(1:end-1),nrs,...
    %                                 lnLambdas(1:end-1),deltas(1:end-1),op.Zs(1:end-1));


    figure(4+ss.figureOffset)
    clf;

    subplot(2,1,1)
    hold on;

    hsText = {'dn_r/dt - Total'};
    hs(1) = plot(sp.timestepsForPlot,dnrdts,'Linestyle','-','Color','k');
    if ss.calculateForOnlyPrimaries
        hsText = [hsText, {'dn_r/dt - Primaries','dn_r/dt - Avalanche'}];
        hs(2) = plot(sp.timestepsForPlot,dnrdtsPrim,'Linestyle','--','Color',colors(2,:));
        hs(3) = plot(sp.timestepsForPlot,dnrdts-dnrdtsPrim,'Linestyle','-.',...
            'Color',colors(1,:));
    end
    hs = [hs, plot(sp.timestepsForPlot,dreicerEstimate+avalancheEstimate,'Linestyle',':',...
        'Color',colors(6,:),'lineWidth',2)];
    hsText = [hsText,'Estimate - Dreicer+Avalanche'];
    if ss.calculateForOnlyPrimaries
        hsText = [hsText, {'Estimate - Dreicer','Estimate - Avalanche'}];
        hs(5) = plot(sp.timestepsForPlot,dreicerEstimate,'Linestyle',':',...
            'Color',colors(5,:),'lineWidth',2);
        hs(6) = plot(sp.timestepsForPlot,avalancheEstimate,'Linestyle',':',...
            'Color',[1,0.8,0],'lineWidth',2);
    end

    ylabel('Runaway generation rate (per reference collision time)');

    sourceTimes = find(sourceRecomputations);
    lims = get(gca,'YLim');
    for i = 2:length(sourceTimes)
        hSource = plot(sp.timesteps(sourceTimes(i))*ones(1,2),lims,'--k','LineWidth',0.5);
        uistack(hSource,'bottom');
        if i == 2
            %Add it to the legend
            hs = [hs hSource];
            hsText = [hsText,'Source recomputations'];
        end
    end

    legend(hs,hsText,'Location','Best');



    subplot(2,1,2)
    hold on;

    plot(sp.timestepsForPlot,nrs,'Linestyle','-','Color','k');
    if ss.calculateForOnlyPrimaries
        plot(sp.timestepsForPlot,nrsPrim,'Linestyle','--','Color',colors(2,:));
        plot(sp.timestepsForPlot,nrs-nrsPrim,'Linestyle','-.','Color',colors(1,:))
    end

    dreicerEstimate = cumsum(sp.dts.*dreicerEstimate);
    avalancheEstimate = cumsum(sp.dts.*avalancheEstimate);
    plot(sp.timestepsForPlot,dreicerEstimate+avalancheEstimate,'Linestyle',':',...
        'Color',colors(6,:),'lineWidth',2);
    if ss.calculateForOnlyPrimaries
        plot(sp.timestepsForPlot,dreicerEstimate,'Linestyle',':',...
            'Color',colors(5,:),'lineWidth',2);
        plot(sp.timestepsForPlot,avalancheEstimate,'Linestyle',':',...
            'Color',[1,0.8,0],'lineWidth',2)
    end

    xlabel('time');% [1/nu_{ee,ref}]');
    ylabel('Runaway density');
    reLims = get(gca,'Ylim');
    %     set(gca,'Ylim',[0,reLims(2)]);
    if ss.calculateForOnlyPrimaries
        legend('n_r/n - Total','n_r/n - Primaries','n_r/n - Avalanche',...
            'Estimate - Dreicer+Avalanche','Estimate - Dreicer',...
            'Estimate - Avalanche','Location','best');
    else
        legend('n_r/n - Total','Estimate - Dreicer+Avalanche','Location','best');
    end

end


if ss.showSourceTerms && ss.initialDistribution ~= 1 %Particle and heat sources
    figure(10+ss.figureOffset)
    clf;

    if sp.noTimeDep
        subplot(1,4,1)
        plot(sp.timesteps,fDensities/fDensities(1));
        xlabel('time');
        ylabel('Density (normalized to the density of a Maxwellian with the same density)');

        subplot(1,4,3)
        plot(sp.timesteps,fHeat/fHeat(1));
        xlabel('time');
        ylabel('Heat (normalized to heat of a Maxwellian at the same temperature)');

    else
        subplot(1,4,1)
        plot(sp.timesteps,fDensities./(sqrt(pi)*sp.nBars/4));
        xlabel('time');
        ylabel('Density (normalized to the density of a Maxwellian with the same density)');

        subplot(1,4,3)
        plot(sp.timesteps,fHeat./(3*sqrt(pi)*sp.nBars.*sp.veBars.^2/8));
        xlabel('time');
        ylabel('Heat (normalized to the heat of a Maxwellian at the same temperature)');
    end

    subplot(1,4,2)
    plot(sp.timesteps,pFactors);
    xlabel('time');
    ylabel('Particle source factor');

    subplot(1,4,4)
    plot(sp.timesteps,hFactors);
    xlabel('time');
    ylabel('Heat source factor');
end


if ss.showCurrentsPlot


    %     if exist('timestepsForPlot','var')
    times = sp.timestepsForPlot;
    %     else
    %         times = op.timesteps;
    %     end
    fracRE = nrs;
    fracRECurrent = currents;

    figure(12+ss.figureOffset)
    clf;

    subplot(3,1,1)
    %     plot(times,nes/(sqrt(pi)/4),...
    %          times,fracRE,'-',...
    %          times,nExternals,'--',...
    %          times,fracRE+nExternals,'-.');
    plot(times,fracRE,'-',...
        times,nExternals,'--',...
        times,fracRE+nExternals,'-.');
    legText = {'RE','External','RE+External'};
    if ss.fastParticleDefinition
        hold on;
        plot(times,nfs,':');
        legText = [legText,'Fast'];
    end
    legend(legText,'Location','best')
    xlabel(sp.timeLabelText);
    ylabel('Runaway fraction');

    subplot(3,1,2)
    plot(times,fracRECurrent,'-',...
        times,externalCurrents,'--',...
        times,fracRECurrent+externalCurrents,'-.');
    if ss.fastParticleDefinition
        hold on;
        plot(times,fastCurrents,':');
    end
    legend(legText,'Location','best')
    xlabel(sp.timeLabelText);
    ylabel('Current fraction');


    subplot(3,1,3)
    plot(times,fracRECurrent.*totCurrentDens,'-',...
        times,externalCurrents.*totCurrentDens,'--',...
        times,(fracRECurrent+externalCurrents).*totCurrentDens,'-.');
    if ss.fastParticleDefinition
        hold on;
        plot(times,fastCurrents.*totCurrentDens,':');
    end
    hold on;
    plot(times,totCurrentDens,'-k');
    legText = [legText,'Total'];

    legend(legText,'Location','best')
    xlabel(sp.timeLabelText);
    ylabel('Current [A/m^2 ]');

end


if ss.showSourceMomentsPlot %Plot properties of the avalanche source
    hFig = figure(14+ss.figureOffset);
    set(hFig,'name','Moments of avalanche source');
    clf;

    if ss.sourceMode == 5
        nRows = 4;
        nCols = 2;
        doSink = 1;
        idDens = 1;
        idMom = 3;
        idEne = 5;
    else
        nRows = 3;
        nCols = 1;
        doSink = 0;
        idDens = 1;
        idMom = 2;
        idEne = 3;
    end

    %Density
    subplot(nRows,nCols,idDens);
    plot(sp.timestepsForPlot,sourceDensity,'-');
    if doSink
        hold on;
        plot(sp.timestepsForPlot,sinkDensity,'--');
        legend('Source','Sink','Location','Best');
        %         legend('Source','Sink','Total runaway source','Location','Best');
    end
    xlabel('Time');
    ylabel('Density');

    %Momentum
    subplot(nRows,nCols,idMom);
    plot(sp.timestepsForPlot,-sourceMomentum,'-');
    if doSink
        hold on;
        plot(sp.timestepsForPlot,-sinkMomentum,'--');
    end
    xlabel('Time');
    ylabel('Momentum');

    %Energy
    subplot(nRows,nCols,idEne);
    plot(sp.timestepsForPlot,sourceEnergy,'-');
    if doSink
        hold on;
        plot(sp.timestepsForPlot,sinkEnergy,'--');
    end
    xlabel('Time');
    ylabel('Energy');

    %Relative values
    if doSink
        subplot(nRows,nCols,2);
        clrs = get(gca,'ColorOrder');
        semilogy(sp.timestepsForPlot,(sourceDensity-sinkDensity)./sourceDensity,'-','Color',clrs(1,:));
        hold on;
        semilogy(sp.timestepsForPlot,-(sourceMomentum-sinkMomentum)./sourceMomentum,'--','Color',clrs(2,:));
        semilogy(sp.timestepsForPlot,(sourceEnergy-sinkEnergy)./sourceEnergy,'-.','Color',clrs(3,:));
        legend('Density','Momentum','Energy','Location','Best');
        semilogy(sp.timestepsForPlot,-(sourceDensity-sinkDensity)./sourceDensity,':','Color',clrs(1,:));
        semilogy(sp.timestepsForPlot,(sourceMomentum-sinkMomentum)./sourceMomentum,':','Color',clrs(2,:));
        semilogy(sp.timestepsForPlot,-(sourceEnergy-sinkEnergy)./sourceEnergy,':','Color',clrs(3,:));
        xlabel('Time');
        ylabel('(Q_{so}-Q_{si})/Q_{so}');
        title('Source conservational properties');


        subplot(nRows,nCols,4)
        semilogy(sp.timestepsForPlot,(reSourceDensity-reSinkDensity)./nrs,'-','Color',clrs(1,:));
        hold on;
        semilogy(sp.timestepsForPlot,-(reSourceMomentum-reSinkMomentum)./reMomentum,'--','Color',clrs(2,:));
        semilogy(sp.timestepsForPlot,(reSourceEnergy-reSinkEnergy)./reEnergy,'-.','Color',clrs(3,:));
        legend('Density (n_S/n_r)','Momentum (p_S/p_r)',...
            'Energy (\gamma_S/\gamma_r)','Location','Best');
        semilogy(sp.timestepsForPlot,-(reSourceDensity-reSinkDensity)./nrs,':','Color',clrs(1,:));
        semilogy(sp.timestepsForPlot,(reSourceMomentum-reSinkMomentum)./reMomentum,':','Color',clrs(2,:));
        semilogy(sp.timestepsForPlot,-(reSourceEnergy-reSinkEnergy)./reEnergy,':','Color',clrs(3,:));
        xlabel('Time');
        title('Normalized total RE source');
    end

end

% The plot showing parameter time dependence is described where inputs are
% processed, as is the distribution evolution plot, the avalanche source
% contour plot, and the grid-in-MeV plot. The plot checking the accuracy of
% grid rescalings is described where initial distributions are interpreted.
% The plot showing the different terms in the collision operator is
% described inside the time loop.


if ss.showSpeciesPlot
    PlotSpecies(ss.species, sp.timeLabelText)
end

%%%%%%%%%%%%%%%%%%%%%%
% Return the result: %
%%%%%%%%%%%%%%%%%%%%%%

%************************************************************************taken from GenerateKnockOnMatrix function


        showPlots = 0;
        if showPlots
            Fplot1 = (p.^2.*pWeights)'*M_source(1:Ny,1:Ny);
            Fplot2 = (gamma' >= (2*gamma_m-1)) .* p'.^3 ./ gamma' .* sigma0' .* pWeights';

            intFunc1 = pWeights'.*Fplot1;
            intFunc1(isnan(intFunc1)) = 0;
            intFunc2 = pWeights'.*Fplot2;
            intFunc2(isnan(intFunc1)) = 0;

            %    Fplot3 = (p.^2.*sqrt(1+p.^2).*pWeights)'*M_source(1:Ny,1:Ny);
            %    Fplot2 = (gamma' >= (2*gamma_m'-1)) .* p'.^3 ./ gamma' .* sigma0Energy' .* pWeights';
            figure(998)
            clf
            subplot(2,1,1)
            plot(p,(Fplot1-Fplot2)./Fplot2,'k')
            axis([0 p(end) -0.1 0.1])
            subplot(2,1,2)
            plot(p,(cumsum(intFunc1) - cumsum(intFunc2))./cumsum(intFunc1),'k')
            set(gca,'xlim',[0 p(end)])



            figure(999)
            clf
            semilogy(p,(p(:).^2.*pWeights(:))'*M_source(1:Ny,1:Ny),'k')
            hold on
            semilogy(p,-2*(p(:).^2.*pWeights(:))'*M_sink1(1:Ny,1:Ny),'--r')
            semilogy(p,-2*(p(:).^2.*pWeights(:))'*M_sink2(1:Ny,1:Ny),':b')
            axis([0.723 0.733 1e-7 1e-3])
        end



        %     fprintf('Calculated knock-on matrix in %0.1f seconds!\n',toc(totalTime))
        pause(1e-6)

        %     % Plot density conservation curve
        %     dm = (pWeights(:).*p(:).^2)'*M(1:Ny,1:Ny);
        %     dm_m = (dm(3:end)+2*dm(2:end-1)+dm(1:end-2))/4;
        %     dm_s = (pWeights(:).*p(:).^2)'*M_source(1:Ny,1:Ny);
        %     dm_sm = (dm_s(3:end)+2*dm_s(2:end-1)+dm_s(1:end-2))/4;
        %     figure(983)
        %     clf
        %     subplot(2,1,1)
        %     plot(p(2:end-1),dm_m./dm_sm,'k')
        %     subplot(2,1,2)
        %     plot(p(1:end),cumsum(dm)/sum(dm_s),'k')

        %     A = whos;
        %     memoryUsages = [A.bytes]' / 1e6;
        %     names = {A.name}';
        %     memIndices = find(memoryUsages > 20);
        %     names(memIndices)
        %     memoryUsages(memIndices)
        %
        %     sprintf('Total memory usage: %0.1f MB',sum(memoryUsages))



%***************************************************************************taken from row 1412 in middle of main for loop
    %Update some plots
    if isUpdateTimestep(iteration)
        if ss.showDistEvolutionPlot
            if sp.noTimeDep
                if ss.calculateForOnlyPrimaries
                    PlotInstantaneousDist(hDistEvoAxes,y,soln,ss.Ny,ss.Nxi,1,...
                        wMinForSource,deltaRef,sp.nBars(2),...
                        sp.veBars3(2),sp.timesteps(iteration),idTail,...
                        ss.yMaxForPlot,ss.FMinForPlot,solnPrim);
                else
                    PlotInstantaneousDist(hDistEvoAxes,y,soln,ss.Ny,ss.Nxi,1,...
                        wMinForSource,deltaRef,sp.nBars(2),...
                        sp.veBars3(2),sp.timesteps(iteration),idTail,...
                        ss.yMaxForPlot,ss.FMinForPlot);
                end
            else
                if ss.calculateForOnlyPrimaries
                    PlotInstantaneousDist(hDistEvoAxes,y,soln,ss.Ny,ss.Nxi,1,...
                        wMinForSource,deltaRef,sp.nBars(iteration),...
                        sp.veBars3(iteration),sp.timesteps(iteration),idTail,...
                        ss.yMaxForPlot,ss.FMinForPlot,solnPrim);
                else
                    PlotInstantaneousDist(hDistEvoAxes,y,soln,ss.Ny,ss.Nxi,1,...
                        wMinForSource,deltaRef,sp.nBars(iteration),...
                        sp.veBars3(iteration),sp.timesteps(iteration),...
                        idTail,ss.yMaxForPlot,ss.FMinForPlot);
                end
            end
        end
        if ss.showAvalancheSourcePlot
            PlotAvalancheSource(hAvaSourceFig,ss,y,sourceVector,sinkVector,...
                sp.timesteps(iteration),wMinForSource,deltaRef);
        end
        if ss.showVarPlot && ~sp.parametersHaveChanged(iteration) %Move the progress bar
            set(hProgressLines(:),'XData',...
                [sp.timesteps(iteration),sp.timesteps(iteration)]);
            drawnow;
        end
        if ss.showCollOpPartPlot
            [fPF,~] = SumLegModesAtXi1(fPPart*fMinus1,ss.Ny,ss.Nxi);
            [tPF,~] = SumLegModesAtXi1(tPPart*fMinus1,ss.Ny,ss.Nxi);

            figure(9+ss.figureOffset)
            clf;
            %            semilogy(x,fPF,x,tPF,x,fPF+tPF);
            plot(x,fPF,x,tPF,x,fPF+tPF);
            xlim([0,5])
            xlabel('x');
            ylabel('F');
            title(['Parts of the collision operator acting on f, shown at v_{perp} = 0 at t = ',num2str(sp.timesteps(iteration)),' reference collision times']);
            legend('C_{ee}^{fp}(f)', 'C_{ee}^{tp}(f)','C_{ee}(f)');
        end
    end

%**************************************************************************taken from row 1584
        if ss.showVarPlot %Update the progress bar/matrix recomputation plot
            if gridWasExtended
                plot(hProgressAx,[sp.timestepsForPlot(iteration),sp.timestepsForPlot(iteration)],...
                    progressAxLim,'-','Color',MyGreen,'linewidth',2);
            else
                plot(hProgressAx,[sp.timestepsForPlot(iteration),sp.timestepsForPlot(iteration)],...
                    progressAxLim,'-k','linewidth',2);
            end

            uistack(hProgressLines(end),'top');

            set(hProgressLines(:),'XData',...
                [sp.timestepsForPlot(iteration),sp.timestepsForPlot(iteration)]);
            drawnow;
        end

%**************************************************************************taken from 1625

        if ss.showCollOpPartPlot
            fPPart = BuildMatrix(ss,sp,iteration,x,y,dxdy,x2,y2,...
                ddy,d2dy2,ddx,d2dx2,deltaRef,'fppart');

            tPPart = BuildMatrix(ss,sp,iteration,x,y,dxdy,x2,y2,...
                ddy,d2dy2,ddx,d2dx2,deltaRef,'tppart');
        end



%******************* plot functions from the end
    function PlotInstantaneousDist(hDistEvo,y,soln,Ny,Nxi,runId,...
            wCrit,deltaRef,nBar,veBar3,time,idTail,...
            yMaxForPlot,FMinForPlot,varargin)

        colors = {'b','r','g','c'};
        colorsP = {[1,0.8,0],'m','b','r'};
        warning off MATLAB:Axes:NegativeDataInLogAxis

        yCrit = wCrit/deltaRef; %Lower boundary of the runaway region

        %Positive and negative here has been reversed since the code was
        %written. Negative X here thus means the direction of the runaway
        %tail.
        cla(hDistEvo);

        [FNegativeX,FPositiveX] = SumLegModesAtXi1(soln,Ny,Nxi);
        xBig = -[fliplr(-y), y]';
        F = [flipud(FNegativeX); FPositiveX];
        semilogy(hDistEvo,xBig, F, 'Color',colors{runId});
        hold(hDistEvo,'on');
        plotHandle = semilogy(hDistEvo,xBig, -F, ':','Color',colors{runId},'LineWidth',2);
        set(get(get(plotHandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        if size(varargin) == 1 %If there is a primary-only-distribtuion
            [FNegativeX,FPositiveX] = SumLegModesAtXi1(varargin{1},Ny,Nxi);
            FPrim = [flipud(FNegativeX); FPositiveX];
            semilogy(hDistEvo,xBig, FPrim,'--','Color',colorsP{runId});
            plotHandle = semilogy(hDistEvo,xBig, -FPrim, ':','Color',colorsP{runId},'LineWidth',2);
            set(get(get(plotHandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end

        semilogy(hDistEvo,xBig, nBar/veBar3 *exp(-xBig.^2/veBar3^(2/3)), '-.r');
        semilogy(hDistEvo,xBig, exp(-xBig.^2), '--k');
        %         semilogy(hDistEvo,xBig, nBar/veBar3 * exp(-xBig.^2'/veBar2) .*
        %         (xBig.^2'/veBar2 - 3/2),'--m'); %This is the heat source term
        xlim(hDistEvo,[-10,yMaxForPlot])
        ylim(hDistEvo,[FMinForPlot, max(max(F),1)])
        if isreal(yCrit) && yCrit > 0
            semilogy(hDistEvo,[yCrit,yCrit],[FMinForPlot, max(max(F),1)],'-k');
            if size(varargin) == 1
                legT = {'Distribution at v_{perp}=0','Distribution - Primaries',...
                    'Instantaneous Maxwellian','Reference Maxwellian',...
                    'Boundary of runaway region'};
            else
                legT = {'Distribution at v_{perp}=0',...
                    'Instantaneous Maxwellian','Reference Maxwellian',...
                    'Boundary of runaway region'};
            end
        else
            if size(varargin) == 1
                legT = {'Distribution at v_{perp}=0','Distribution - Primaries',...
                    'Instantaneous Maxwellian','Reference Maxwellian'};
            else
                legT = {'Distribution at v_{perp}=0',...
                    'Instantaneous Maxwellian','Reference Maxwellian'};
            end
        end
        if idTail
            semilogy(hDistEvo,[y(idTail),y(idTail)],[FMinForPlot, max(max(F),1)],'--m');
            legT = [legT, 'Boundary of fast particle region'];
        end
        legend(hDistEvo,legT);
        title(hDistEvo,['Total distribution for v_{perp} = 0 at t = ',num2str(time),' reference collision times']);
        xlabel(hDistEvo,'\gamma v_{||} / v_{th,ref}')
        ylabel(hDistEvo,'F')
        %         pause(0.2)
        drawnow;

    end

    function PlotAvalancheSource(hAvaSourceFig,ss,yGrid,sourceVector,...
            sinkVector,time,wCrit,deltaRef)

        doTreatSink = any(sinkVector);

        nPoints = 80;
        aspect = 0.6;
        fractionOnNegSide = 0.2;

        %Set up a grid
        paraGrid = linspace(0,ss.yMax,nPoints);
        perpGrid = paraGrid(1:round(aspect*nPoints));
        negSide = -fliplr(paraGrid(2:round(fractionOnNegSide*nPoints)));
        paraGrid = [negSide,paraGrid];
        [par2D,perp2D] = meshgrid(paraGrid,perpGrid);

        source2D = GetDistributionPointYY(par2D,perp2D,yGrid,sourceVector,ss.Nxi,ss.Ny);
        if doTreatSink
            sink2D = GetDistributionPointYY(par2D,perp2D,yGrid,sinkVector,ss.Nxi,ss.Ny);
        end

        %Determine sensible contours to use
        mx = log10(max(source2D(:)));
        sourceContours = linspace(mx-10,mx,30);

        %Prepare yCrit and yM circles
        yCrit = wCrit/deltaRef;
        yM = ss.yCutSource;

        nPCircle = 30;
        angles = linspace(0,pi,nPCircle);
        circlePara = cos(angles);
        circlePerp = sin(angles);

        %Plot it
        clf;
        if doTreatSink
            hSink = subplot(2,1,2,'Parent',hAvaSourceFig);

            contourf(hSink,par2D,perp2D,log10(sink2D),sourceContours,'EdgeColor','none');
            colorbar(hSink)
            hold(hSink,'on');
            hCrit = plot(hSink,yCrit*circlePara,yCrit*circlePerp,'--k');
            hM = plot(hSink,yM*circlePara,yM*circlePerp,'-.r');
            legend(hSink,[hCrit,hM],{'Boundary of RE region','Source cutoff'})
            xlabel(hSink,'y_{||}')
            ylabel(hSink,'y_\perp');
            title(hSink,['Avalanche Sink term (log_{10}(-S)) at t = ',num2str(time)]);
            axis(hSink,'equal');

            hSou = subplot(2,1,1,'Parent',hAvaSourceFig);
        else
            hSou = axes('Parent',hAvaSourceFig);
        end

        contourf(hSou,par2D,perp2D,log10(source2D),sourceContours,'EdgeColor','none');
        colorbar(hSou)
        hold(hSou,'on');
        hCrit = plot(hSou,yCrit*circlePara,yCrit*circlePerp,'--k');
        %         hM = plot(hSou,yM*circlePara,yM*circlePerp,'-.r');
        %         legend([hCrit,hM],'Boundary of RE region','Source cutoff')
        legend(hSou,hCrit,'Boundary of RE region')
        xlabel(hSou,'y_{||}')
        ylabel(hSou,'y_\perp');
        title(hSou,['Avalanche Source term (log_{10}(S)) at t = ',num2str(time)]);
        axis(hSou,'equal');

        drawnow;
    end

    function [vals,ts] = ExpandForPlot(vals,ts,allTimes)
        if isscalar(vals) && isscalar(ts)
            vals = vals*[1,1];
            ts = allTimes([1,end]);
        end
    end
%********************************* taken from functions at the end, used for runaway plots
    function dnrdt = EstimateDreicerGeneration(EOverED,Z,nueeBars)
        %Estimate the Dreicer runaway generation with the formula used in
        %GO (see Fehér, PPCF, 2011). The growth rate is per reference collision time
        EDE = 1./EOverED;
        exp1 = 3*(Z+1)/16;
        exp2 = -0.25*EDE-sqrt((Z+1).*EDE);

        dnrdt = 0.25*3*sqrt(pi)*nueeBars .* EDE.^exp1 .* exp(exp2);
        %Note that nr in CODE is actually nr/n, so for a comparison the
        %normalization cancels the n in this formula.

    end

    function dnrdt = EstimateAvalancheGeneration(EOverEc,nr,lnLambda,delta,...
            Z,nueeBars)
        %Estimate the avalanche runaway generation with the formula used in
        %GO (see Fehér, PPCF, 2011). The growth rate is per reference collision time
        preFactor = 0.25*3*sqrt(pi)*nr.*delta.^3.*(EOverEc-1).*nueeBars./lnLambda;
        rootArg = pi./(3*(Z+5));
        lastTerm = 4*pi*(Z+1).^2 ./ ( 3*(Z+5).*(EOverEc.^2+3) );

        dnrdt = preFactor.*sqrt( rootArg./(1-1./EOverEc+lastTerm) );
    end

%******************************* taken from functions at the end, used for plotting avalance source

    function fOut = GetDistributionPointYY(yPara,yPerp,yGrid,f,Nxi,Ny)
        % Calculates the value of F. The distribution is returned for points in
        % momentum space given by the corresponding arrays yPara & yPerp.
        % The value of f is calculated using F(y,xi) = sum_l F_l(y) P_l(xi).

        ys = sqrt(yPara.*yPara + yPerp.*yPerp);
        xi = yPara./ys;

        %Flatten the ys and xi arrays (make them row vectors)
        outSize = size(ys);
        ys = ys(:)';
        xi = xi(:)';

        %Generate a matrix of values of the Legendre polynomials at the given xi.
        %legendres has the structure
        % [ leg(1,xi(1)), leg(1,xi(2)), ...,  leg(1,xi(nXi)) ]
        % [ leg(2,xi(1)), leg(2,xi(2)), ...,  leg(2,xi(nXi)) ]
        % [     ...     ,      ...    , ...,        ...      ]
        % [leg(Nxi,xi(1)),leg(Nxi,xi(2)),...,leg(Nxi,xi(nXi))]

        %Using the purpose-built routine (this is ~100 times faster)
        legendres = LegendrePolynomials(Nxi-1,xi);
        fLMatrix = reshape(f,Ny,Nxi);

        %Interpolate all the modes to the given y in one call. fLAtY has the
        %structure:
        % [ f0(y(1)),  f1(y(1)),..., fNxi(y(1)) ]
        % [ f0(y(2)),  f1(y(2)),..., fNxi(y(2)) ]
        % [   ...   ,    ...   ,...,     ...    ]
        % [ f0(y(nY)),f1(y(nY)),..., fNxi(y(nY))]
        fLAtY = interp1(yGrid',fLMatrix,ys);

        %Calculate each mode at the given points.
        fLs = legendres' .* fLAtY;

        %Sum over the modes
        plusMin = ones(1,Nxi);
        plusMin(2:2:end) = -1;
        fOut = sum(fLs*diag(plusMin),2);

        fOut = reshape(fOut,outSize(1),outSize(2));

        %Ignore negative values. A cutoff can also be introduced
        cutoff = 0;%1e-15;
        fOut(fOut < cutoff) = 0;%cutoff;

    end


  %******************************* taken from functions at the end, used for plotting species

        function PlotSpecies(species, timelabel) %fix labels etc later...
        figure(15+ss.figureOffset)
        clf
        for iSpecies = 1:size(species.nj,1)
            plot(species.times, species.nj(iSpecies,:))
            hold on
        end
        legend(strcat('Z=', num2str(species.ZAtomicNumber), ' Z0 = ',  num2str(species.Z0NetCharge)))
        title('Time evolution of species densities')
        xlabel(timelabel)
        ylabel('n')
        axSpecies = gca;
        axSpecies.YLim = [0 axSpecies.YLim(2)*1.1];
    end

