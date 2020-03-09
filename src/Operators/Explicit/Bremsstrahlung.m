classdef Bremsstrahlung < ExplicitOperator
    % BOLTZMANN generates the bremstrahlung matrix and then returns the operatorF
    %WITHOUT the sourcePreFactor in original CODE (which is to account for
    %different solution methods (and thus should be done in the solver))
    % basically does the things GenerateMatricesForBoltzmannOperators did
    % minus the knock on matrix plus doing the line:
    % boltzmannMatrix = nBar*(1+this.plasma.Z(iteration))*bremsMatrix

    properties
            NxiUsed = 0     
            pUsed = 0;
            NyInterpUsed = 0
            useScreeningUsed = -1
            bremsModeUsed = -1;
            speciesUsed %todo smart comparision with only at current timestep
            
            
            nBarUsed = 0
            ZUsed = 0
            deltaUsed = 0
            lnLambdaRefUsed = 0
    end
    
    methods
        function this = Bremsstrahlung(state,eqSettings)
            this@ExplicitOperator(state,eqSettings)
        end

        function matrixHasChanged = generateOperatorMatrix(this,iteration)
            matrixHasChanged = 1;
            
            nBar = this.state.physicalParams.nBars(iteration);
            Z = this.state.physicalParams.Z(iteration);
            delta = this.state.physicalParams.deltas(iteration);
            lnLambda = this.state.reference.lnLambdaRef;

            sz = this.state.momentumGrid.matrixSize;
            if this.NxiUsed == this.state.momentumGrid.Nxi && ...
                isequal(this.pUsed, this.state.reference.deltaRef*this.state.momentumGrid.y(:)) && ...
                this.NyInterpUsed  == this.eqSettings.NyInterp && ...
                this.useScreeningUsed == this.eqSettings.useScreening && ...
                this.bremsModeUsed  == this.eqSettings.bremsMode && ...
                ~((this.eqSettings.bremsMode == 4)&&(this.eqSettings.useScreening == 2)&&~isequal(this.speciesUsed, this.state.physicalParams.species)) %only checks species equality if bremsMode 4 and source 2(only time species is used)
                %everything is same (except for some constants maybe
                if this.nBarUsed == nBar &&...
                        this.ZUsed == Z &&...
                        this.deltaUsed == delta &&...
                        this.lnLambdaRefUsed == lnLambda
                    convFac = (delta/this.deltaUsed)^3/(lnLambda/this.lnLambdaRefUsed)*(nBar/this.nBarUsed)*(1+Z)/(1+this.ZUsed);
                    this.operatorMatrix = convFac*this.operatorMatrix;
                else
                    matrixHasChanged = 0;
                    return
                end
                
            else
                if this.eqSettings.bremsMode >=2
                    switch this.eqSettings.bremsMode
                        case 2
                            smallKOperator= 1;
                            angularDepMode = 'full';
                        case 3
                            smallKOperator= 1;
                            angularDepMode = 'delta';
                        case 4
                            smallKOperator= 0;
                            angularDepMode = 'delta';
                        case 5
                            smallKOperator= 0;
                            angularDepMode = 'delta';
                            %Also adds energy diffusion to Fokker-Planck operator
                        otherwise
                            error('Something went wrong -- invalid bremsMode.');
                    end
                    if this.eqSettings.useScreening && (this.eqSettings.bremsMode~= 4 || this.eqSettings.useScreening~= 2|| numel(this.state.physicalParams.species.times)>1)
                        warning('Screened bremsstrahlung only works with bremsMode 4, useScreening= 2 and no time-dependence in species')
                    end
                    alpha = 1/137.036;
                    preFac = 3/(16*sqrt(pi)) * alpha*delta^3/lnLambda;
                    [brMat, ~, ~] = this.GenerateBremsstrahlungMatrix(angularDepMode,smallKOperator);
                    brMat = preFac * brMat;

                else
                    brMat = sparse(sz,sz);
                end
                this.operatorMatrix = nBar*(1+Z)*brMat; %bremsMatrix effectively includes ntot/nfree already
            end
            
            this.NxiUsed            = this.state.momentumGrid.Nxi;
            this.pUsed              = this.state.reference.deltaRef*this.state.momentumGrid.y(:);
            this.NyInterpUsed       = this.eqSettings.NyInterp;
            this.useScreeningUsed   = this.eqSettings.useScreening;
            this.bremsModeUsed      = this.eqSettings.bremsMode;
            this.speciesUsed        = this.state.physicalParams.species;
            
            
            this.nBarUsed           = this.state.physicalParams.nBars(iteration);
            this.ZUsed              = this.state.physicalParams.Z(iteration);
            this.deltaUsed          = this.state.physicalParams.deltas(iteration);
            this.lnLambdaRefUsed    = this.state.reference.lnLambdaRef;
        end
        
    end

    methods
        function brMat = GenerateBremsMatrix(this,iteration)
            delta = this.state.physicalParams.deltas(iteration);
            lnLambda = this.state.reference.lnLambdaRef;

            sz = this.state.momentumGrid.matrixSize;

            %Bremsstrahlung operator
            if this.eqSettings.bremsMode >=2
                switch this.eqSettings.bremsMode
                    case 2
                        smallKOperator= 1;
                        angularDepMode = 'full';
                    case 3
                        smallKOperator= 1;
                        angularDepMode = 'delta';
                    case 4
                        smallKOperator= 0;
                        angularDepMode = 'delta';
                    case 5
                        smallKOperator= 0;
                        angularDepMode = 'delta';
                        %Also adds energy diffusion to Fokker-Planck operator
                    otherwise
                        error('Something went wrong -- invalid bremsMode.');
                end
                if this.eqSettings.useScreening && (this.eqSettings.bremsMode~= 4 || this.eqSettings.useScreening~= 2|| numel(this.state.physicalParams.species.times)>1)
                    warning('Screened bremsstrahlung only works with bremsMode 4, useScreening= 2 and no time-dependence in species')
                end
                alpha = 1/137.036;
                preFac = 3/(16*sqrt(pi)) * alpha*delta^3/lnLambda;
                [this.M, this.M_source, this.M_sink] = this.GenerateBremsstrahlungMatrix(angularDepMode,smallKOperator);
                brMat = preFac * [this.M this.M_source, this.M_sink];

            else
                brMat = sparse(sz,sz);
            end
        end

        function [M,M_source,M_sink] = GenerateBremsstrahlungMatrix(this,mode,smallK)
            yWeights = this.state.momentumGrid.yWeights;
            Ny       = this.state.momentumGrid.Ny;
            Nxi      = this.state.momentumGrid.Nxi;
            p        = this.state.reference.deltaRef*this.state.momentumGrid.y(:);
            pAll =     this.state.reference.deltaRef*this.state.momentumGrid.y(:);
            pMax     = p(end);
            gamma    = sqrt(1+p.^2);
            pMin_Ind = 2; % accounts only for those collisions where the outgoing electron
            % has energy >= the first non-zero grid value. (as diff-cross
            % sec is singular at photon energy = theoretical maximum, although
            % the limit is well-defined, and zero. This is equivalent.)

            kMinOld  = gamma/1000;

            upperI      = find((kMinOld+gamma)>gamma(end),1)-1;%Ny-1;
            pOld        = p;
            %gammaOld    = gamma;
            p           = p(pMin_Ind:upperI);
            gamma       = sqrt(1+p.^2);
            Np          = length(p);

            kMin     = gamma/1000; %yielding an approximate 0.1% error in the energy loss.

            pMin        = p(1);
            %gammaMin    = gamma(1);
            % pkMin_Ind   = find(gammaOld >= (gammaMin + kMin),1);
            % kMin        = gammaOld(pkMin_Ind) - gammaMin;
            % fprintf(' kMin = %0.2d, kMin/EMax = %0.2d \n pMin = %0.2d, pMin/pMax = %0.2d \n',...
            %                                            kMin,kMin/(gamma(end)-1),pMin,pMin/p(end))

            P           = repmat(p,1,this.eqSettings.NyInterp);
            Gamma       = sqrt(1+P.^2);
            P1          = pMax*ones(size(P));
            WS1         = zeros(size(P));

            p_int_min   = sqrt((gamma+kMin).^2-1);
            parfor K = 1:Np                     % Interpolating to a Gauss-Legendre-
                % quadrature grid (using the lgwt func)
                [ps,ws]  = lgwt(this.eqSettings.NyInterp,p_int_min(K),pMax);
                P1(K,:)  = fliplr(ps);
                WS1(K,:) = fliplr(ws);
            end

            Gamma1 = sqrt(1 + P1.^2);

            M0 = P1.^3 .*WS1 ./ (P.*Gamma.*Gamma1);
            clear('WS1')
            if strcmp(mode,'delta')
                if this.eqSettings.useScreening == 2 && this.eqSettings.bremsMode == 4
                    M0 = M0.* this.d1sigma_BetheHeitler_screened(P1,Gamma1-Gamma,this.state.physicalParams.species);
                else
                    M0 = M0.* this.d1sigma_BetheHeitler(P1,Gamma1-Gamma);
                end
                clear('Gamma1','Gamma','P')

                t_memSave = tic;
                blockSize = 64;

                Lmat     = repmat(0:(Nxi-1),Ny,1);
                oneMat   = ones(Ny,Nxi);
                Jmat     = repmat((1:Ny)',1,Nxi) + Ny*Lmat;
                J0       = Jmat(:);
                M_source = spalloc(Ny*Nxi,Ny*Nxi,round(Np^2*Nxi/2));

                maxBlockNumber = floor(Np/blockSize);
                remainderSize = mod(Np,blockSize);

                for blockNumber = 1:maxBlockNumber
                    I   = cell(blockSize,1);
                    S   = cell(blockSize,1);

                    P_interp_array  = cell(blockSize,1);
                    M0_array       = cell(blockSize,1);
                    p_int_mins      = p_int_min((1:blockSize) + blockSize*(blockNumber-1));
                    for R = 1:blockSize
                        K = R + blockSize*(blockNumber-1);
                        P_interp_array{R} = P1(K,:);
                        M0_array{R}      = M0(K,:);
                    end


                    parfor R = 1:blockSize
                        %                 sprintf('R=%d, ',R)
                        K = R + blockSize*(blockNumber-1);
                        Imat    = K*oneMat + Ny*Lmat;
                        I{R}    = Imat(:) + pMin_Ind-1;
                        S{R}    = zeros(size(I{R}));
                        if pMax > p_int_mins(R)
                            interpMatrix    = this.BuildInterpolationMatrix(pOld,P_interp_array{R});
                            %                         ML0t            = repmat(squeeze(M0_array{R}),Nxi,1);
                            MLt             = (repmat(squeeze(M0_array{R}),Nxi,1) * interpMatrix)';
                            S{R} = MLt(:);
                        end
                    end
                    I = cell2mat(I);
                    S = cell2mat(S);
                    J = repmat(J0,blockSize,1);
                    M_source = M_source + sparse(I,J,S,Ny*Nxi,Ny*Nxi);
                end

                if remainderSize > 0
                    I   = cell(remainderSize,1);
                    S   = cell(remainderSize,1);
                    for R = 1:remainderSize
                        K = R + blockSize*maxBlockNumber;
                        Imat    = K*oneMat + Ny*Lmat;
                        I{R}    = Imat(:)+pMin_Ind-1;
                        S{R}    = zeros(size(I{R}));
                        if pMax > p_int_min(K)
                            interpMatrix    = this.BuildInterpolationMatrix(pOld,P1(K,:));
                            MLt             = (repmat(squeeze(M0(K,:)),Nxi,1) * interpMatrix)';
                            S{R} = MLt(:);
                        end
                    end
                    I = cell2mat(I);
                    S = cell2mat(S);
                    J = repmat(J0,remainderSize,1);
                    M_source = M_source + sparse(I,J,S,Ny*Nxi,Ny*Nxi);
                end
                toc(t_memSave)
                clear('S,I,J,MLt,PL,ML0t,P1')


            elseif strcmp(mode,'full')
                I   = zeros(Np*Ny*Nxi,1);  % The full I,J will represent Nxi solid
                J   = I;
                Lmat = repmat(0:(Nxi-1),Ny,1);
                Jmat = repmat((1:Ny)',1,Nxi) + Ny*Lmat;
                for K = 1:Np
                    % This is indeed a pretty funky way of indexing the matrix.
                    % I am building it row-by-row for all Legendre modes simultaneously.
                    indices     = (K-1)*Ny*Nxi + (1:(Ny*Nxi));
                    Imat        = K*ones(Ny,Nxi) + Ny*Lmat;
                    I(indices)  = Imat(:) + pMin_Ind-1;
                    J(indices)  = Jmat(:);
                end
                S = cell(Ny,1);
                MLs = zeros(Nxi,Np,this.eqSettings.NyInterp);
                for L = 0:Nxi-1
                    numMintervals = 5;
                    numNintervals = 5;
                    Mranges = [0:2,ceil(this.eqSettings.NyInterp/(10*numMintervals)),round((1:numMintervals)*this.eqSettings.NyInterp/numMintervals)];%5,20,40,80,140,numInterpPoints];
                    Nranges = [0:5,ceil(Np/(10*numNintervals)),round((1:numNintervals)*Np/numNintervals)];%200,800,1100,Ny];
                    NN = length(Nranges);
                    MM = length(Mranges);

                    reltol = 5e-4;
                    abstol = 1e-6;
                    for N = 1:NN-1
                        for M = 1:MM-1
                            Ns = (1+Nranges(N)):Nranges(N+1);
                            Ms = (1+Mranges(M)):Mranges(M+1);
                            integrand = @(theta) sin(theta).*this.LastLegendrePolynomial(L,cos(theta)) ...
                                .*this.diffCrossSec(P(Ns,Ms),P1(Ns,Ms),theta);
                            if length(Ns)==1 && length(Ms)==1
                                MLs(L+1,Ns,Ms) = integral(@(theta)integrand(theta),0,pi,'reltol',reltol,'abstol',abstol);
                            else
                                MLs(L+1,Ns,Ms) = integral(@(theta)integrand(theta),0,pi,'ArrayValued',true,'reltol',reltol,'abstol',abstol);
                            end
                        end
                    end
                    MLs(L+1,:,:) = M0.*squeeze(MLs(L+1,:,:));
                end
                MLs(isnan(MLs)) = 0;
                parfor K = 1:Np
                    interpMatrix = this.BuildInterpolationMatrix(pOld,P1(K,:));
                    MLt     = squeeze(MLs(:,K,:)) * interpMatrix;
                    MLt     = MLt';
                    S{K}    = MLt(:);
                end

                S = cell2mat(S);
                clear('MLs')
                M_source = sparse(I,J,S,Ny*Nxi,Ny*Nxi);
                clear('I','J','S')

            end
            pause(1e-6)

            % CONSIDER INSERTING BACK THIS COMMENTED-OUT PIECE. Current method shown to be problematic for knock-on case.
            % totalCrossSec = zeros(size(p));
            % for K = pkMin_Ind:Ny
            %     totalCrossSec(K) = integral(@(k)this.d1sigma_BetheHeitler(pOld(K),k),...
            %                                 kMin,gammaOld(K)-gammaMin);
            % end
            % sinkFunc = -pOld./gammaOld .* totalCrossSec;
            % M_sink   = spdiags(repmat(sinkFunc(:),Nxi,1),0,Ny*Nxi,Ny*Nxi);
            % This could possibly be improved by choosing the quadrature weights as a
            % simpson quadrature
            sinkFunc = -((pAll(:).^2.*yWeights(:))'*M_source(1:Ny,1:Ny))./(pAll(:).^2.*yWeights(:))';
            sinkFunc(isnan(sinkFunc)) = 0;
            sinkFunc(isinf(sinkFunc)) = 0;
            M_sink   = spdiags(repmat(sinkFunc(:),Nxi,1),0,Ny*Nxi,Ny*Nxi);
            M        = M_source + M_sink;

            if smallK
                % NEEDS TO BE FIXED: use proper k_min involving n_e instead of 7e-10
                Msmall = this.GeneratePitchOperator(pOld,Nxi,7e-10,kMinOld);
                M      = M + Msmall;
            end

            outParams.pMin = pMin;
            outParams.kMin = kMinOld;
        end

    end
    %************************* helpfunctions for building the matricies ********************************************
    methods
        function d1sigma= d1sigma_BetheHeitler(~,p,k)
            E = sqrt(1+p.^2);
            E1 = E-k;
            p1 = sqrt(E1.^2-1);

            PreFactor = p1./(p.*k);

            E0E = E.*E1;
            p0p = p.*p1;

            Eps  = 2*log(E1+p1);
            Eps0 = 2*log(E+p);
            L = 2*log( (E0E + p0p - 1)./k );

            Term1 = 4/3 - 2*E0E.* (p1.^2 + p.^2)./(p0p.^2) ...
                + Eps0.*E1./p.^3 + Eps.*E./p1.^3 - Eps.*Eps0./(p0p);
            Term2 = L.*( 8*E0E./(3*p0p) + k.^2./(p0p.^3) .*( E0E.^2 + p0p.^2 ) );
            Term3 = L.* k./(2*p0p) .*( (E0E + p.^2).*Eps0./p.^3 ...
                - (E0E+p1.^2).*Eps./p1.^3 + 2*k.*E0E./p0p.^2 );

            d1sigma = PreFactor .* ( Term1 + Term2 + Term3 );
            d1sigma(k<=0) = 0;
            d1sigma(k>=E-1) = 0;
        end

        function d1sigma_screened = d1sigma_BetheHeitler_screened(this,p0,k,species)
            % FormFactor defined with q =  2sin(theta/2)*(p/mc)

            E0 = sqrt(p0.^2+1);
            E  = E0-k;
            p  = sqrt(E.^2-1);

            g1 = 1 + (E./E0).^2;
            g2 = -2*E./(3*E0);

            q0 = p0-p-k;
            q = @(x) q0 + (1-q0).*x; % variable substitution to be able to use 'arrayvalued'.

            iterationTMP = 1;
            ZMinusFSquared = this.calculateAveragedFormFactorBrems(species, iterationTMP);

            nj_t = species.nj(:,iterationTMP); %density of each species
            Zeffne = (species.Z0NetCharge.^2)' * nj_t;
            Zeff_fullne = (species.ZAtomicNumber.^2)' * nj_t;

            FI1 = Zeff_fullne +   (1-q0).*integral(@(x)  ZMinusFSquared(q(x)) .* (q(x)-q0).^2./q(x).^3 , 0,1,'arrayvalued',true);
            FI2 = 5/6*Zeff_fullne + (1-q0).*integral(@(x)  ZMinusFSquared(q(x)) .* (q(x).^3 + 3*q(x).*q0.^2.*( 1 - 2*log(q(x)./q0)) -4*q0.^3)./q(x).^4,0,1,'arrayvalued',true);

            d1sigma_screened = 4./k.*(g1.*FI1 + g2.*FI2)/Zeffne;%n_efree*Zeff: normalize to this
            % to cancel a factor of Zeff somewhere else in CODE...

        end
        function C = GeneratePitchOperator(~,p,Nxi,kMin,kMax)
            Nxi0 = Nxi;
            if Nxi>201
                warning('Pitch operator only defined up to L = 200. Setting C_L = 0 for L>200.')
                Nxi  = 200;
            end
            if p(end) > 10000 %|| (p(1) < 1)
                warning('Pitch operator only defined in the range p = [1,10000]mc. Setting C = 0 for p>10000mc and p<mc.')
            end

            mod = log(kMax./kMin) / log(100); % used kMax/kMin = 100 to generate data

            Ny      = length(p);
            CLMat   = load('matfiles/CLMat3.mat');
            CLmat   = CLMat.CL;
            P1Mat   = load('matfiles/P12.mat');
            p1      = P1Mat.p1;

            % We interpolate in p*C on a x-log scale, as this is a more slowly varying
            % function.
            p1Mat   = spdiags(p1(:),0,length(p1),length(p1));
            % pMat    = spdiags(p(:),0,Ny,Ny);
            % modMat  = spdiags(mod(:),0,Ny,Ny);

            mod = mod(:);
            p = p(:);
            v = p./sqrt(1+p.^2);

            preFactorVec = - mod.*v./p;
            preFactorVec(1) = 0;
            preFactor = spdiags(preFactorVec,0,Ny,Ny);

            CL      = preFactor * interp1(log10(p1),(CLmat(1:Nxi-1,:)*p1Mat)',log10(p),'pchip',0);

            % filling with zeros L=0 or L > 200
            ZerosForLeq0 = zeros(Ny,1);
            ZerosForLgtr200 = zeros(Ny*(Nxi0-Nxi),1);

            Cdiag   = [ZerosForLeq0; CL(:); ZerosForLgtr200];
            C       = spdiags(Cdiag,0,Ny*Nxi0,Ny*Nxi0);
        end

        function [ ZMinusFSquared ] = calculateAveragedFormFactorBrems(this, species, iteration)
            % q is in units of p/mc - so add a factor of alpha! (but no factor of 2)
            % and F is normalized to Z
            % uses useScreening = 2

            n = 3/2;
            ALPHA = 1/137.036;
            ZMinusFSquared = @(q) 0;

            nj_t = species.nj(:,iteration); %density of each species
            Z0NetCharge = species.Z0NetCharge;

            load('DFTdata/aTFDFT', 'aTFDFT');
            atomicSymbolList = {'He','Be','C','N','Ne','Ar','Xe','W'};
            atomicNumberList = [ 2,   4,   6,  7,  10,  18,  54,  74];


            for ispecies = 1:numel(Z0NetCharge)
                Z = species.ZAtomicNumber(ispecies);
                Z0 = Z0NetCharge(ispecies);
                Ne = Z-Z0;
                nj = nj_t(ispecies);

                if Ne>0 % otherwise F_j = 0
                    %Translate Z,Z0 into NnX (eg. Ar1) and check if data exists
                    if any(Z==atomicNumberList)
                        atomicSymbol = atomicSymbolList{Z==atomicNumberList};
                    else
                        atomicSymbol = 'NaN';
                    end
                    speciesName = [atomicSymbol num2str(Z0)];
                    if (~isfield(aTFDFT, speciesName))
                        error('no data on chosen species.')
                    end

                    a_lengthscale = aTFDFT.(speciesName);
                    ZMinusFSquared = @(q) ZMinusFSquared(q) + nj*(Z-Ne./(1+(q*a_lengthscale/ALPHA).^n)).^2;
                end
            end
        end


        function W = diffCrossSec(this,p,p1,theta)
            % Evaluates dsigma/dkdOmega

            costhetas = cos(theta);
            gamma1  = sqrt(1 + p1.^2);
            gamma   = sqrt(1 + p.^2);
            k       = gamma1 - gamma;

            Gamma1Gamma = gamma1.*gamma;
            Gamma1GammaSq = Gamma1Gamma.^2;
            kSq  = k.^2;
            p1Sq = p1.^2;
            p1Cu = p1Sq.*p1;
            pSq  = p.^2;
            pCu  = pSq.*p;
            p1p  = p1.*p;
            p1pSq = p1p.^2;
            p1pQa = p1pSq.^2;
            Gamma1Sq = gamma1.^2;
            GammaSq  = gamma.^2;

            l       = log(p+gamma);
            l1      = log(p1+gamma1);
            lambda  = Gamma1Gamma - p1p.*costhetas - 1;
            lambdaSq = lambda.^2;
            lambdaCu = lambdaSq.*lambda;
            lambdaQa = lambdaSq.^2;
            SqrtLambda = sqrt(lambdaSq + 2*lambda);
            Term1 = (2*Gamma1Gamma + (Gamma1Sq + GammaSq - 1).*lambda - lambdaSq) ...
                ./(kSq.*lambdaSq.*SqrtLambda).*log(1+lambda+SqrtLambda);
            Term2   = - (2*Gamma1Gamma - lambda)./(kSq.*lambdaSq);
            Term3   = - 3*(Gamma1GammaSq - 1).^2./(lambdaSq.*p1pQa);
            Term4   = (4*(Gamma1GammaSq - Gamma1Gamma + 1) ...
                - Gamma1Gamma.*(p1Sq + pSq) + (Gamma1Sq + GammaSq + Gamma1Gamma - 1).*lambda) ...
                ./(2*lambdaSq.*p1pSq);
            Term5   = (2*(Gamma1Gamma - 1)./lambdaCu - kSq./lambdaQa) ...
                .*(2*(Gamma1Sq + GammaSq - Gamma1Gamma).*p1pSq + 3*kSq.*(gamma1 + gamma).^2) ...
                ./p1pQa;

            Term6   = l./pCu .* gamma.*(1+2*GammaSq)./(lambdaSq.*pSq);
            Term7   = l./pCu.*(2*GammaSq.^2 + 2*p1pSq + gamma1.*(gamma1 + gamma) ...
                - (Gamma1Gamma + pSq).*lambda)./(2*k.*lambdaSq);
            Term8   = l./pCu.*gamma.*(2*(Gamma1Gamma - 1)./lambdaCu - kSq./lambdaQa) ...
                .*(2*gamma1.*pSq - 3*k.*GammaSq)./(k.*pSq);
            Term9   = l1./p1Cu.*gamma1.*(1+2*Gamma1Sq)./(lambdaSq.*p1Sq);
            Term10  = - l1./p1Cu.*(2*Gamma1Sq.^2 + 2*p1pSq + gamma.*(gamma1 + gamma) ...
                - (Gamma1Gamma + p1Sq).*lambda)./(2*k.*lambdaSq);
            Term11  = -l1./p1Cu.*gamma1.*(2*(Gamma1Gamma - 1)./lambdaCu - kSq./lambdaQa) ...
                .*(2*gamma.*p1Sq + 3*k.*Gamma1Sq)./(k.*p1Sq);



            W = 2*p.*k./p1.*(Term1 + Term2 + Term3 + Term4 + Term5 ...
                + (Term6 + Term7 + Term8) ...
                + (Term9 + Term10 + Term11) );

            W(p>p1) = 0; % Does not satisfy energy conservation

        end

        function outPl = LastLegendrePolynomial(this,l,x)
            PLs = LegendrePolynomials(l,x(:)');
            outPl = PLs(end,:);
            outPl = reshape(outPl,size(x));
        end

        function K = BuildInterpolationMatrix(this,x,xp)
            Nxp = length(xp);
            Nx  = length(x);
            I = 1:Nxp;
            J = zeros(size(I));
            S = ones(size(I));

            iterator = 2;
            J(1) = 1;
            for K = 1:Nx
                if x(K) > xp(1)
                    iterator = K-1;
                    J(1)     = iterator;
                    break
                end
            end
            for K = 2:Nxp
                while  x(iterator) < xp(K)
                    iterator = iterator+1;
                end
                J(K) = iterator-1;
            end

            D    = sparse([I,1],[J,Nx-1],[S,0]);
            Dbar = sparse([I,1],[J,Nx-1],[S.*xp,0]);

            c = 1./(x(2:end)-x(1:end-1));
            lenc = length(c);
            c1 = x(2:end).*c;
            c2 = -x(1:end-1).*c;
            A = spdiags([-c(:),c(:)],[0,1],lenc,lenc+1);
            B = spdiags([c1(:),c2(:)],[0,1],lenc,lenc+1);
            K = (Dbar*A + D*B);
        end
    end
end
