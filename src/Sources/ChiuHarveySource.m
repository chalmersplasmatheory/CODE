classdef ChiuHarveySource < Source
    %ChiuHarvey like source
    properties
    end

    methods
        function this = ChiuHarveySource(state,eqSettings)
            this@Source(state,eqSettings)
        end
    end

    methods
        function sourceVector = getSourceVec(this,f,iteration)
            
            y = this.state.momentumGrid.y;
            y2 = this.state.momentumGrid.y2;
            x = this.state.momentumGrid.x;
            Nxi = this.state.momentumGrid.Nxi;
            Ny = this.state.momentumGrid.Ny;
            deltaRef = this.state.reference.deltaRef;
            nueeRef = this.state.reference.nueeRef;
            nRef = this.state.reference.nRef;
            f0 = f(1:Ny);
            EOverEc = this.state.physicalParams.EOverEc(iteration);
            nBar = this.state.physicalParams.nBars(iteration);
            w = y * deltaRef;
            matrixSize = Nxi*Ny;
            sourceMode = this.eqSettings.sourceMode;
            neTotalOverneFree = this.state.physicalParams.neTotalOverneFree(iteration);
            yWeights = this.state.momentumGrid.yWeights;
            
            if abs(EOverEc) > 1
                gammaC = sqrt(abs(EOverEc)/(abs(EOverEc)-1) );
                wMinForSource = sqrt(4*gammaC*(gammaC-1)); % only count runaways with
                                                     % kinetic energy twice the
                                                     % critical runaway energy
            else
                wMinForSource = inf;
                gammaC = inf;
            end

            mask = w >= wMinForSource;

            %%%%%%%%%prel for chiuHarvey:
            r0SqG = (2.817940e-13)^2; %Gaussian
            cG = 2.997925e+10; %Gaussian
            nRefG = 1e-6*nRef; %Gaussian
            preFactorConstant = pi*r0SqG*cG*nRefG/nueeRef;
            cHSourcePreFactor = (preFactorConstant*2*deltaRef^3 * x./y2)';
            cHSourcePreFactor(1) = 0; %Remove NaN, this point is not used anyway (no runaways with 0 velocity...)
            gamma = y./x;
            gamma(1) = 1;
            gMax = gamma(end);

            %Initialize an empty vector of the right size
            sourceVector = zeros(matrixSize,1);


            % y-integral implementation

            %%% Build chSourceMatrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %We can restrict the size of cHSourceMatrix, since not all
            %combinations of (gamma,gamma'_1) are allowed.
            if sourceMode == 3
                idYMaxForCHSource = find((gamma-1)>(0.5*(gMax-1)),1)-1;
                %From the condition on the maximum energy of the outgoing electron
                %for the maximum possible incoming electron (yMax)

                % Build Moller cross section ----------------------------------
                [gammaIn2D,gamma2D] = meshgrid(gamma,gamma(1:idYMaxForCHSource));
                idsToDiscard = (gammaIn2D-1) < 2*(gamma2D-1);
                mollerCS = CalculateMollerCS(gamma2D,gammaIn2D);
                mollerCS(idsToDiscard) = 0;
            end

            %!!! This will not work with very large grids! We need to do
            %    some kind of sparse treatment. !!!

            % -------------------------------------------------------------

            % Build P_l matrix ----------------------------------
            XiFunc = @(g,g1) sqrt( (g-1).*(g1+1)./((g+1).*(g1-1)) );
            xi = XiFunc(gamma2D,gammaIn2D);
            xi(idsToDiscard) = 0;
            xiSize = size(xi); %We need it as a vector to calculate the Legendre polynomials
            xiVec = -xi(:)';
            %The minus sign is to reverse the direction of the source, due
            %to the unusual definition of xi in CODE.

            legPols  = LegendrePolynomials(Nxi-1,xiVec);
            %Transform into the form P_l(gamma,gamma'_1,L). Note that mode
            %L is accessed using P_l(:,:,L+1)!
            legPols = reshape(legPols',xiSize(1),xiSize(2),Nxi);
            % -------------------------------------------------------------

            %Combine into an array with (gamma,gamma'_1,L)
            cHSourceMatrix = legPols .* repmat(mollerCS,[1,1,Nxi]);
            %     %This takes up too much memory for large problems, but doing it at
            %     %runtime for each mode takes too long. We will have to think of something.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            cHMomentumVector = (y.^3)./gamma;
            

            %Chiu-Harvey-like source
            idYFirstAboveCrit = find(mask,1);
            wMinCH = 2 * sqrt(gammaC*(gammaC-1));
            intMask = w >= wMinCH;
            idYMinCHInt = find(intMask,1);
            cHIntVector = neTotalOverneFree*yWeights.*cHMomentumVector.*intMask;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Chiu-Harvey-like source, y-int implementation
            intOnGrid = zeros(Ny,1);
            idYForSource = idYFirstAboveCrit:idYMaxForCHSource;

            for L = 0:Nxi-1
                 indices = L*Ny+(1:Ny);

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
            
            this.sourceVector = sourceVector;
        end
    end
end

