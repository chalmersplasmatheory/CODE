classdef CollisionOperator < ImplicitOperator
    %CODESWITCHABLEOPERATOR builds CollisionOperator term (not electric field term)
    %Basically does everything
    %done under the switch variable CollisionalOperator, generate matrix is
    %the same as BuildMatrix function except removed Efield Term

    
    %Bremsmode and SourceMode is left with this.*.* in code since it seems
    %they should not belong here (but might rightfully do so) and is therefore
    %not written in the same standards to raise questions
    properties
        deltaUsed = 0
        veBarUsed = 0
        nueeBarUsed = 0
        lnLambdaUsed = 0
        ZUsed = 0
        NxiUsed = 0
        RelLimitUsed = 0
        yUsed = 0
        xUsed = 0
        collisionOperatorUsed = -1 %to make sure we rebuilt
        useInelasticUsed = 0
        useEnergyDependentLnLambdaScreeningUsed = 0
        sourceModeUsed = 0
        neTotalOverneFreeUsed = 0
        yMaxBCUsed = 0
        artificialDissipationWidthUsed = 0
        artificialDissipationStrengthUsed = 0
        useScreeningUsed = 0
        speciesUsed = 0
        yCutSourceUsed = 0
    end
    
    methods
        function this = CollisionOperator(state,eqSettings)
            %CollisionOperator Construct an instance of this class
            %   Detailed explanation goes here
            this@ImplicitOperator(state,eqSettings);
        end

        function matrixHasChanged = generateOperatorMatrix(this,iteration)
            % initiate values
            
            delta       = this.state.physicalParams.deltas(iteration);
            veBar       = this.state.physicalParams.veBars(iteration);
            veBar3      = veBar^3;
            nueeBar     = this.state.physicalParams.nueeBars(iteration);
            lnLambda    = this.state.physicalParams.lnLambdas(iteration);
            Z           = this.state.physicalParams.Z(iteration);
            species     = this.state.physicalParams.species;
            neTotalOverneFree = this.state.physicalParams.neTotalOverneFree(iteration);
            
            Ny = this.state.momentumGrid.Ny;
            Nxi = this.state.momentumGrid.Nxi;
            yMax = this.state.momentumGrid.yMax;
            sqrtpi = sqrt(pi);

            collisionOperator = this.eqSettings.collisionOperator;
            useInelastic = this.eqSettings.useInelastic;
            useEnergyDependentLnLambdaScreening = this.eqSettings.useEnergyDependentLnLambdaScreening;
            useScreening = this.eqSettings.useScreening;
            
            switch collisionOperator
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
                y2 = this.state.momentumGrid.x2;
                ddy = this.state.momentumGrid.ddx;
                d2dy2 = this.state.momentumGrid.d2dx2;
                dxdy = ones(size(this.state.momentumGrid.dxdy));
                deltaRef = 0;
            else
               y = this.state.momentumGrid.y;
               y2 = this.state.momentumGrid.y2;
               ddy = this.state.momentumGrid.ddy;
               d2dy2 = this.state.momentumGrid.d2dy2;
               dxdy = this.state.momentumGrid.dxdy;
               deltaRef = this.state.reference.deltaRef;
            end

            w = deltaRef*y;
            x = this.state.momentumGrid.x;
            x2 = this.state.momentumGrid.x2;
            gamma = this.state.momentumGrid.gamma;
            ddx = this.state.momentumGrid.ddx;
            d2dx2 = this.state.momentumGrid.d2dx2;


            veBar2 = veBar*veBar;
            expx2 = exp(-x2/veBar2);
            
            %check if rebuild is necessary
            if  abs((this.deltaUsed - delta)/delta)<1e-13 && ...
                    this.veBarUsed == veBar && ...
                    this.nueeBarUsed == nueeBar && ...
                    this.lnLambdaUsed == lnLambda && ...
                    this.ZUsed == Z && ...
                    this.NxiUsed == Nxi && ...
                    this.RelLimitUsed == takeNonRelLimit && ...
                    isequal(this.yUsed, y) && ...
                    isequal(this.xUsed, x) && ...
                    this.collisionOperatorUsed == collisionOperator && ...
                    this.useInelasticUsed == useInelastic && ...
                    this.useEnergyDependentLnLambdaScreeningUsed == useEnergyDependentLnLambdaScreening && ...
                    this.sourceModeUsed == this.eqSettings.sourceMode && ...
                    this.neTotalOverneFreeUsed == neTotalOverneFree && ...
                    this.yMaxBCUsed == this.state.momentumGrid.yMaxBoundaryCondition && ...
                    this.artificialDissipationWidthUsed == this.state.momentumGrid.artificialDissipationWidth && ...
                    this.artificialDissipationStrengthUsed == this.state.momentumGrid.artificialDissipationStrength && ...
                    this.useScreeningUsed == useScreening && ...
                    isequal(this.speciesUsed, species) && ...
                    this.yCutSourceUsed == this.eqSettings.yCutSource
                matrixHasChanged = 0;
                return
            end
            
            disp('CO builds');
            %Build mateix

            if collisionOperator == 0
                Theta = delta^2/2;
                Psi0exp = w.*integral(@(s) exp(-(sqrt(1+w.^2.*s.^2)-1)/Theta)...
                     ./sqrt(1+w.^2.*s.^2),0,1,'arrayvalued',true);
                Psi1exp = w.*integral(@(s) exp(-(sqrt(1+w.^2.*s.^2)-1)/Theta),0,1,'arrayvalued',true);

                if Theta > 1e3/511e3
                    K2exp   = besselk(2,1/Theta).*exp(1/Theta);
                else
                    K2exp   = sqrt(pi*Theta/2).*(1 + 15*Theta/8 + 15*7*Theta.^2/128);
                end


                Psi     = Theta*(gamma.^2.*Psi1exp - Theta*Psi0exp + ...
                          (Theta*gamma-1).*w.*exp(-(gamma-1)/Theta))./(w.^2.*K2exp);


                PsiPrime =  Theta*(2*w.*Psi1exp + gamma.^2.* exp(-(sqrt(1+w.^2)-1)/Theta) ...
                            - Theta*exp(-(sqrt(1+w.^2)-1)/Theta)./sqrt(1+w.^2) ...
                            + (Theta*(w.^2+gamma.^2) - gamma.^3 + w.^2/Theta )...
                            .*exp(-(gamma-1)/Theta)./gamma)./(w.^2.*K2exp) - 2./w .*Psi;
                PsiPrime = PsiPrime*deltaRef./dxdy;
                %%%%% Note that there is a factor deltaRef in the below "missing"
                %%%%% (or well, make sure that you get it correctly in the final
                %%%%% expression)
                nuD_thermal = ( (w.^2.*gamma.^2+Theta^2).*Psi0exp + ...
                    Theta*(2*w.^4-1).*Psi1exp + gamma.*Theta.*(1+Theta*(2*w.^2-1)).*...
                    w.*exp(-(gamma-1)/Theta)).*x.*deltaRef./(gamma.*w.^3*K2exp);


            else
                erfs = erf(x/veBar);
                % Psi is the Chandrasekhar function:
                Psi = (erfs-2/(veBar*sqrtpi)*x.*expx2) * veBar2 ./ (2*x2);
                PsiPrime = 2/(veBar*sqrtpi)*expx2 - 2*Psi./x; %Mine
                nuD_thermal = erfs - Psi + 0.5*veBar2*(deltaRef^4)*x2;
            end
            NL = 2;


            % NOTE: This version of CODE uses a significantly different system of indexing than previous versions.
            % In particular, every Legendre mode of F has all Ny values of x stored.

            % Order of rows in the matrix and right-hand side:
            % --------------------
            % for L = 0:(Nxi-1)
            %   Impose boundary condition at y=0
            %   for iy = 2:(Ny-1)
            %     Impose kinetic equation
            %   Impose boundary condition at yMax
            % Enforce total density = n
            % Enforce total pressure = n*T

            % Order of columns in the matrix, corresponding to rows in the solution vector:
            % --------------------
            % for L = 0:(Nxi-1)
            %   for iy = 1:Ny
            %     Value of F
            % Particle source
            % Heat source

            % Predict roughly how many nonzero elements will be in the sparse
            % matrix. This speeds up the code by eliminating the need to reallocate
            % memory during matrix construction:
            this.predictedNNZ = 2 * Nxi * nnz(abs(ddy))+ Nxi*(nnz(abs(ddy)) + nnz(abs(d2dy2))) ... % Test particle part of collision operator
                + max([NL,0])*Ny*Ny ... % Field term in Fokker-Planck operator %It was min in Matt's code. Must be an error?
                + 4*Ny; % Constraints and sources

            if collisionOperator == 3
                xPotentials = y';
                ddxPotentials = ddy;
                d2dx2Potentials = d2dy2;
                NxPotentials = Ny;
                xMax = y(end);
                expxfp2 = exp(-y2/veBar2);
                xfp2 = y2;
            else
                xPotentials = x';
                ddxPotentials = ddx;
                d2dx2Potentials = d2dx2;
                NxPotentials = Ny;
                xMax = x(end);
                xfp2 = x2;
                expxfp2 = expx2;
            end

            x2 = x.*x;

            xWith0s = [0; xPotentials(2:(end-1)); 0];
            M21 = 4*pi*spdiags(xWith0s.^2,0,Ny,Ny);
            M32 = -2*spdiags(xWith0s.^2,0,Ny,Ny);
            LaplacianTimesX2WithoutL = spdiags(xPotentials.^2,0,Ny,Ny)*d2dx2Potentials + 2*spdiags(xPotentials,0,Ny,Ny)*ddxPotentials;

            % Momentum part of test particle term
            if (~useInelastic  && ~useEnergyDependentLnLambdaScreening && ~(this.eqSettings.sourceMode==5))
                xPartOfCECD = 3*sqrt(pi)/4 * (spdiags((veBar3*nueeBar.*Psi./x)',0,Ny,Ny)*d2dy2   + ...
                    spdiags((veBar3*nueeBar.*((dxdy./x).* PsiPrime + 2*Psi./(x.*y) ...
                    - (Psi./x2).*dxdy + 2/veBar2*Psi))',0,Ny,Ny)*ddy + ...
                    veBar*spdiags((nueeBar.*(2*dxdy.*PsiPrime + 4*Psi./y))',0,Ny,Ny));
            else % change nuS due to inelastic collisions, or change the Coulomb logarithm.
                p = deltaRef*y;

                %calculate energy dependence of Coulomb logarithm
                if useEnergyDependentLnLambdaScreening == 1 || this.eqSettings.sourceMode==5
                    [  lnLambdaeeHat,dlnLambdaeeHatdp,~,boltzCorrectionlnLambdaeeHat, boltzCorrectiondlnLambdaeeHatdp] = ...
                        this.CalculateEnergyDependentlnLambda( p, y, lnLambda);

                else
                    lnLambdaeeHat = ones(size(y));
                    dlnLambdaeeHatdp = zeros(size(y));
                    boltzCorrectionlnLambdaeeHat = zeros(size(y));
                    boltzCorrectiondlnLambdaeeHatdp = zeros(size(y));
                end
                %calculate effect of inelastic collisions
                if useInelastic
                    [ Ghat, dGhatdp ] = this.CalculateInelasticEnhancement(species, ...
                        p, y, lnLambda, useInelastic, iteration );
                else
                    Ghat = zeros(size(y));
                    dGhatdp = zeros(size(y));
                end

                %terms that are not affected by this
                unchangedTerms  = 3*sqrt(pi)/4 * ...
                    (spdiags((veBar3*nueeBar.*Psi./x)',0,Ny,Ny)*d2dy2   + ...
                    spdiags((veBar3*nueeBar.*((dxdy./x).* PsiPrime + 2*Psi./(x.*y) ...
                    - (Psi./x2).*dxdy))',0,Ny,Ny)*ddy);

                % with source 5, correction to avoid double counting inelastic (bound) and elastic (free) electron nu_s
                lnLambdaeeHat_nuS    = (lnLambdaeeHat-neTotalOverneFree*boltzCorrectionlnLambdaeeHat);
                dlnLambdaeeHatdp_nuS = (dlnLambdaeeHatdp-neTotalOverneFree*boltzCorrectiondlnLambdaeeHatdp);

                % \nu_S terms: change due to energy-dependent lnLambda and/or
                % effect of inelastic collisions

                nuSTerms = 3*sqrt(pi)/4 * ( ...
                    spdiags((veBar*nueeBar.*( 2*Psi.*lnLambdaeeHat_nuS + Ghat./x2 ))',0,Ny,Ny)*ddy + ...
                    spdiags((veBar*nueeBar.*( 2*dxdy.*PsiPrime.*lnLambdaeeHat_nuS+...
                    2*deltaRef.*dlnLambdaeeHatdp_nuS.*Psi+ 4*Psi.*lnLambdaeeHat_nuS./y+...
                    2*deltaRef.^2./y.*Ghat+ deltaRef./x2.*dGhatdp))',0,Ny,Ny));

                xPartOfCECD =  unchangedTerms+nuSTerms;

            end

            if this.eqSettings.bremsMode == 1
                ln2g = log(2*gamma);
                xPartOfCECD = xPartOfCECD + ...
                    nueeBar*veBar3*deltaRef^2/lnLambda * 3*sqrt(pi)/4 *(1+Z)/(137*pi) ...
                    *( spdiags(((gamma-1).*(ln2g-1/3))',0,Ny,Ny)*ddy ...
                    + spdiags( (2*(gamma-1).*(ln2g-1/3)./y ...
                    + deltaRef^2*x.*(ln2g +2/3 - 1./gamma))',0,Ny,Ny) );
            end

            switch collisionOperator
                case {0,4}
                    % Relativistic test-particle term
                    %Nothing to do
                case {1,2,3}
                    % Full linearized Fokker-Planck operator
                    xPartOfCECD = xPartOfCECD + 3*spdiags((nueeBar.*expxfp2)',0,Ny,Ny);
                    M12IncludingX0 = 3/(2*pi*veBar2)*spdiags((nueeBar.*expxfp2)',0,Ny,Ny);
                    M13IncludingX0 = -3/(2*pi*veBar*veBar3)* spdiags((nueeBar.*xfp2.*expxfp2)',0,Ny,Ny) * d2dx2Potentials;
            end

            % Smooth the distribution functions at yMax:
            if this.state.momentumGrid.yMaxBoundaryCondition == 4
                fakeViscosity = exp((y-yMax)/this.state.momentumGrid.artificialDissipationWidth);
                xPartOfCECD = xPartOfCECD + this.state.momentumGrid.artificialDissipationStrength*spdiags(fakeViscosity',0,Ny,Ny)*d2dy2;
            end
            
            if this.state.momentumGrid.yMaxBoundaryCondition == 4
                dyEnd = y(end)-y(end-1);
                fakeViscosity = exp( (y-yMax)/(this.state.momentumGrid.artificialDissipationWidth*dyEnd) );
                xPartOfCECD = xPartOfCECD + this.state.momentumGrid.artificialDissipationStrength*spdiags(fakeViscosity',0,Ny,Ny)*d2dy2;
            end

            % Initialize arrays for building the sparse matrix:
            this.resetSparseCreator()

            % calculate the (for Boltzmann) L-dependent enhancement of the collision
            % frequency, gives a matrix of size (Lmax,lenght(p)) %NB: starts at
            % L=1!
            if (useScreening || useEnergyDependentLnLambdaScreening || (this.eqSettings.sourceMode == 5) )
                %but before we do all that, we need to change the coulomb
                %logarithm; it is difficult to do afterwards!

                p = deltaRef*y;

                if useEnergyDependentLnLambdaScreening || (this.eqSettings.sourceMode == 5) %
                    % dependent coulomb logarithm
                    [ lnLambdaeeHat,~,lnLambdaeiHat, boltzCorrectionlnLambdaeeHat, ~] = ...
                        this.CalculateEnergyDependentlnLambda( p, y, lnLambda, delta);
                else
                    [lnLambdaeeHat,lnLambdaeiHat] = deal(ones(size(y)));
                    boltzCorrectionlnLambdaeeHat = zeros(size(y));
                end
                % make lnLambdaeiHat a matrix
                Lmax = Nxi-1;
                lnLambdaeiHat = ones(Lmax,1)*lnLambdaeiHat;

                % calculate the screening factor but add the zeroth Legendre mode
                if useScreening
                    screeningFactor = [ones(size(y));
                        lnLambdaeiHat.*this.CalculateScreeningFactor(species,...
                        p,Lmax,lnLambda*lnLambdaeiHat,useScreening, iteration)];
                else
                    screeningFactor = ones(Nxi,length(y));
                end
                % with source 5, correction to avoid double counting elastic (free) electron nu_D
                lnLambdaeeHat_nuD = (lnLambdaeeHat-boltzCorrectionlnLambdaeeHat);
                
            else % do not consider screening
                screeningFactor = ones(Nxi,length(y));
                lnLambdaeeHat_nuD = 1;
            end


            for L=0:(Nxi-1)

                % Energy- and L-dependent deflection frequency:
                nuD = nueeBar*veBar3*3*sqrt(pi)/4.*(Z*screeningFactor(L+1,:) +...
                    lnLambdaeeHat_nuD.*nuD_thermal) ./ (x.*y2);
                % now the Z part has a factor of lnL^ei(p) and the ee part has a
                % prefactor of lnL^ee(p)

                % Add collision operator
                M11 = - (-0.5*spdiags(nuD',0,Ny,Ny)*L*(L+1) + xPartOfCECD);


                if L <= (NL-1) && ~(collisionOperator== 0 || collisionOperator == 4)
                    % Add Rosenbluth potential terms

                    M13 = M13IncludingX0;
                    M12 = M12IncludingX0;

                    M22 = LaplacianTimesX2WithoutL-L*(L+1)*speye(NxPotentials);

                    % Add Dirichlet or Neumann boundary condition for
                    % potentials at x=0:
                    if L==0
                        M22(1,:)=ddxPotentials(1,:);
                    else
                        M22(1,:) = 0;
                        M22(1,1) = 1;
                        M12(:,1) = 0;
                        M13(:,1) = 0;
                    end
                    M33 = M22;

                    % Add Robin boundary condition for potentials at x=xMax:
                    M22(NxPotentials,:) = xMax*ddxPotentials(NxPotentials,:);
                    M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1;

                    M33(NxPotentials,:) = xMax*xMax*d2dx2Potentials(NxPotentials,:) + (2*L+1)*xMax*ddxPotentials(NxPotentials,:);
                    M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1);
                    if L~=0
                        M22(NxPotentials,1)=0;
                        M33(NxPotentials,1)=0;
                    end

                    collisionOperatorMatrix = - (M11 -  (M12 - M13 * (M33 \ M32)) * (M22 \ M21));

                else
                    collisionOperatorMatrix = - M11;
                end

                rowIndices = L*Ny + (2:(Ny-1));
                columnIndices = L*Ny + (1:Ny);
                this.addSparseBlock(rowIndices, columnIndices, collisionOperatorMatrix(2:(Ny-1), 1:Ny))
            end

            this.operatorMatrix = -this.createSparse();
            this.clearSparseResidue()
            
            
            this.deltaUsed = delta;
            this.veBarUsed = veBar;
            this.nueeBarUsed = nueeBar;
            this.lnLambdaUsed = lnLambda;
            this.ZUsed = Z;
            this.NxiUsed = Nxi;
            this.RelLimitUsed = takeNonRelLimit;
            this.yUsed = y;
            this.xUsed = x;
            this.collisionOperatorUsed = collisionOperator;
            this.useInelasticUsed = useInelastic;
            this.useEnergyDependentLnLambdaScreeningUsed = useEnergyDependentLnLambdaScreening;
            this.sourceModeUsed = this.eqSettings.sourceMode;
            this.neTotalOverneFreeUsed = neTotalOverneFree;
            this.yMaxBCUsed = this.state.momentumGrid.yMaxBoundaryCondition;
            this.artificialDissipationWidthUsed = this.state.momentumGrid.artificialDissipationWidth;
            this.artificialDissipationStrengthUsed = this.state.momentumGrid.artificialDissipationStrength;
            this.useScreeningUsed = useScreening;
            this.speciesUsed = species;
            this.yCutSourceUsed = this.eqSettings.yCutSource;
            matrixHasChanged = 1;
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%% methods for screening and other
    %%%%%%%%%%%%%%%%%%%%%%%%%%% functionalities

    methods
       function [ lnLambdaeeHat,dlnLambdaeeHatdp,lnLambdaeiHat, boltzCorrectionlnLambdaeeHat, boltzCorrectiondlnLambdaeeHatdp] = ...
                CalculateEnergyDependentlnLambda(this, p, y, lnLambda0, delta)
            %calculateEnergyDependentlnLambdlnLambda calculates the enhancement of the Coulomb
            %logarithm replacing the termal speed with the relativistic momentum at
            %high energy. Basically, nuS -> nuS*(lnLambdaHat + Ghat). For details,
            % see Linnea's notes.
            %   The formula is a simple interpolation between the energy dependent
            %   expression and the thermal speed expression at low energy. References:
            %   Wesson p. 792 and (lnL0 = lnL^ee) and Solodev-Betti (2008) for
            %   energy dependence
            % OUTPUT: lnLambdaeehat, the derivative dlnLambdaeehatdp and lnLambdaeihat
            % have the same size as y.
            n = 5; % sharpness of transition between 'thermal' form and high energy form

            %first electron-electron coulomb logarithm
            gamma = sqrt(p.^2+1); % because overwritten by some Boltzmann function

            lnLambdaee = lnLambda0*ones(size(p));
            lnLambdaei = lnLambda0*ones(size(p));
            lnLambdaeePrime = 0;
            if this.eqSettings.useEnergyDependentLnLambdaScreening
                r = sqrt(2*(gamma-1))/delta;
                lnLambdaee = lnLambdaee + 1/n *log(1+r.^n);
                lnLambdaeePrime = lnLambdaeePrime  + p./(2*gamma.*(gamma-1)) ./ (1 + r.^(-n) );
                lnLambdaei = lnLambdaei + 1/n *log(1+(2*y).^n);
            end
            if this.eqSettings.sourceMode == 5 % compensates for Boltzmann test-particle operator
                idGammaM = find(y>this.eqSettings.yCutSource,1);
                gammaM = gamma(idGammaM);
                boltzCorrectionlnLambdaee 	   =   (gamma>gammaM).*log(sqrt((gamma-1)/(gammaM-1)));
                boltzCorrectionlnLambdaeePrime =   (gamma>gammaM).*p./(2*gamma.*(gamma-1));
            else
                boltzCorrectionlnLambdaee 	   =   0;
                boltzCorrectionlnLambdaeePrime =   0;
            end

            lnLambdaeeHat = lnLambdaee./lnLambda0;
            lnLambdaeiHat = lnLambdaei./lnLambda0;
            boltzCorrectionlnLambdaeeHat = boltzCorrectionlnLambdaee./lnLambda0;

            dlnLambdaeeHatdp = lnLambdaeePrime./lnLambda0;
            dlnLambdaeeHatdp(1) = 0;
            boltzCorrectiondlnLambdaeeHatdp = boltzCorrectionlnLambdaeePrime./lnLambda0;

            %electron-ion coulomb logarithm
       end

        function [ Ghat, dGhatdp ] = CalculateInelasticEnhancement(this, species, p, y,lnLambda0, useInelastic, iteration )
            %CALCULATEINELASTICENHANCEMENT calculates the enhancement of the
            %electron-electron slowing down frequency (\nu_S^{ee}) due to inelastic
            %collisions with bound electrons around the ion. Basically, nuS ->
            % nuS*(lnLambdaHat + Ghat). For details see Linnea's notes.
            %   INPUTS:
            %       - species: struct with fields nj, Z0NetCharge, ZAtomicNumber,
            %         times:
            %           - nj(jspecies, times) - matrix containing number density of
            %             each species as a function of time (m^{-3})
            %           - ZAtomicNumber - atomic number of each species. row vector,
            %             not a function of time
            %           - Z0NetCharge - net charge for each species. row vector, not a
            %             function of time
            %           - times - time steps. If constant, just set i to 1.
            %       - p = normalized momentum gamma*m*v/(m*c) = delta*y
            %       - LnLambda the coulomb logarithm
            %       - useInelastic: which model to use
            %           - 0: Ignore inelastic collisions
            %           - 1: Bethe formula with values of mean excitation energy
            %               from Sauer et al.
            %	    - 2: RP Rule of thumb (cutoff at y=10 to conserve particles)
            %  OUTPUT: Ghat and the derivative dGhatdp have the same size as y.
            nj_t = species.nj(:,iteration); %density of each species
            delta = this.state.physicalParams.deltas(iteration);
            Z0NetCharge = species.Z0NetCharge;
            ZAtomicNumber = species.ZAtomicNumber;
            Ne = ZAtomicNumber-Z0NetCharge; %number of bound elecrons
            ne_free = Z0NetCharge' * nj_t; %ne_free
            gamma = sqrt(p.^2+1); % Gamma is already a global variable! But now it doesn't run so maybe safest like this?
            n = 5;% sharpness of transition between 'thermal' form and high energy form
            beta = p./gamma;

            zeroVec = zeros(size(p));
            Ghat = zeroVec;
            dGhatdp = zeroVec;
            switch useInelastic
                case 0 %ignore inelastic collisions
                    % stay 0
                case 1 % bethe with Ihat from sauer et al.
                    % load data
                    load('DFTdata/Ihat', 'Ihat');


                    for i = 1:numel(nj_t)%for every species
                        % the contribution to screening for each species:

                        if Ne(i) ==0
                            g = zeroVec;
                            gPrime =zeroVec;

                        else
                            %First, identify the name of this species - translate input into NnX
                            %(eg. Ar1)
                            atomicSymbolList = {'He','Be','C','N','Ne','Ar','Xe','W'};
                            atomicNumberList = [ 2,   4,   6,  7,  10,  18,  54,  74];
                            if any(ZAtomicNumber(i)==atomicNumberList)
                                atomicSymbol = atomicSymbolList{ZAtomicNumber(i)==atomicNumberList};
                            else
                                atomicSymbol = 'NaN';
                            end
                            speciesName = [atomicSymbol num2str(Z0NetCharge(i))];

                            %now calculate g_i!
                            if isfield(Ihat, speciesName)
                                z = p.*sqrt(gamma-1)./(Ihat.(speciesName)); % auxilary variable
                                g = 1/n.*log(1+z.^n)-beta.^2;
                                gPrime = (3*gamma+1)./(2*gamma.*p)./(1+z.^(-n))-...
                                    2*(p./gamma.^4);
                            else
                                %data missing
                                warning('no data on chosen species. Switching to useInelastic = 0 for this species.\n')
                                g = 0;
                                gPrime =0;
                            end
                        end

                        Ghat = Ghat + ...
                            1/(ne_free) * Ne(i) * nj_t(i) .* (g/lnLambda0);
                        dGhatdp = dGhatdp + ...
                            1/(ne_free) * Ne(i) * nj_t(i) .* (gPrime/lnLambda0);


                    end
                case 2
                    ne_bound = (ZAtomicNumber - Z0NetCharge)' * nj_t; %ne_bound
                    GhatHalfBound = 1/2*ne_bound/ne_free;

                    %density conservation + no change to nuD and nuParallel -> use RP rule outside of bulk only
                    RP_CUTOFF = 10; % apply it from here
                    XSHARPNESS = 1; % sharpness of transition
                    GhatMask = 0.5*(1+tanh(XSHARPNESS*(y-RP_CUTOFF)));
                    GhatMaskPrime = 0.5*XSHARPNESS*(sech(XSHARPNESS*(y-RP_CUTOFF))).^2;

                    GhatWithoutLnLambda = GhatHalfBound*GhatMask;
                    dGhatdpWithoutLnLambda = GhatHalfBound*GhatMaskPrime;
                    if useEnergyDependentLnLambdaScreening == 1
                        [ lnLambdaeeHat,dlnLambdaeeHatdp,~ ] = ...
                            this.CalculateEnergyDependentlnLambda( p, y, lnLambda0, delta);
                        Ghat = GhatWithoutLnLambda.*lnLambdaeeHat;
                        dGhatdp = dGhatdpWithoutLnLambda.*dlnLambdaeeHatdp;
                    elseif useEnergyDependentLnLambdaScreening  % two benchmark versions of coulog
                        warning('did not bother implementing this combination of energy-dependend coulomb logarithm and inelastic collisions \n');
                    end

                    %dGhatdp stays zero
                otherwise
                    warning('This switch for useInelastic is not implemented yet. Must be an integer 0-1.\n');
            end

        end

        function [ g ] = CalculateSpeciesScreeningFactor(this, useScreening,Z,Ne,Z0,p,Lmax )
            %CalculateSpeciesScreeningFactor calculates the function g(p): the \nu =
            %(1+g/lnLambda), for each legendre mode. p must be a row vector.
            % output: g, with size (Lmax, Ny)
            %   INPUTS:
            %     species has the fields
            %     nj, Z0NetCharge, ZAtomicNumber, times:
            %     nj(jspecies, times) - matrix containing number density of each
            %       species as a function of time (m^{-3})
            %     ZAtomicNumber - atomic number of each species. row vector, not a
            %     function of time
            %     Z0NetCharge - net charge for each species. row vector, not a
            %     function of time
            %     times - time steps. If constant, just set i to 1.
            %     - p = normalized momentum gamma*m*v/(m*c) = delta*y
            %     - LnLambda the coulomb logarithm. It is now a
            %       (Nxi-1, Ny ) matrix! (p-dependent)
            %     - useScreening: which model to use
            %       - 0: Ignore screening (complete screening)
            %       - 1: Fokker--Planck collision operator, TF model (works for all species)
            %       - 2: Fokker--Planck collision operator, TF-DFT model (works for some
            %			ions)
            %       - 3: Fokker--Planck collision operator, full DFT model (works for some
            %			ions)
            %       - 4: Boltzmann collision operator, TF model
            %       - 5: Boltzmann collision operator, TF-DFT model
            %       - 6: Boltzmann collision operator, full DFT model
            %       - 7: no screening = full penetration
            %           ADD WARNING IF P IS TOO BIG!!!

            n = 3/2;
            ALPHA = 1/137.036;
            pToq = 2/ALPHA; % convert to q from p/mc
            calcg_TFDFT_y = @ (y) (Z^2-Z0^2)/n*log(1+y.^n)-(Z-Z0)^2/n*y.^n./(1+y.^n);
            calcg_TFDFT_p = @(p,a)calcg_TFDFT_y(p*a*pToq);
            calcF_TFDFT_y = @(y) Ne./(1+(y).^n);
            calcF_TFDFT_p = @(p,a)calcF_TFDFT_y(p*a*pToq);

            a_TF = @(Z) 1 * ((9*pi^2)/(128*Z))^(1/3); %a0=1 in normalized units
            aTF_Krillov = @(Z,Ne)((Ne)^2/Z^2*2/pi)^(1/3)*a_TF(Z);

            if Ne==0
                g = zeros(Lmax, length(p)); %no need to do all this
                return
            end

            %__________________________________________________________________________
            %First, find the fieldname - translate input into NnX (eg. Ar1)
            %First replace Z with atomic number.
            atomicSymbolList = {'He','Be','C','N','Ne','Ar','Xe','W'};
            atomicNumberList = [ 2,   4,   6,  7,  10,  18,  54,  74];

            if any(Z==atomicNumberList)
                atomicSymbol = atomicSymbolList{Z == atomicNumberList};
            else
                atomicSymbol = 'NaN';
            end
            ionization = Z-Ne;
            speciesName = [atomicSymbol num2str(ionization)];
            %__________________________________________________________________________

            %load data
            switch useScreening
                case 2 %FP, TFDFT
                    load('DFTdata/aTFDFT', 'aTFDFT');
                case 3 %FP, DFT
                    load('DFTdata/gDFT','gDFT');
                case 5 %Boltz, TFDFT
                    load('DFTdata/aTFDFT', 'aTFDFT');
                    load('DFTdata/gBoltzTFDFT','gBoltzTFDFT');
                case 6 %Boltz, DFT
                    load('DFTdata/FDFT', 'FDFT');
                    load('DFTdata/gBoltzDFT','gBoltzDFT');
            end
            %__________________________________________________________________________

            %If TFDFT or DFT: test if g data exists, otherwise use TF

            if (useScreening == 2 && ~isfield(aTFDFT, speciesName)) ||...
                    (useScreening == 3 && ~isfield(gDFT, speciesName))
                %data missing from FP, DFT or TFDFT, change to FP with Thomas--Fermi
                warning('no data on chosen species. Using Fokker--Planck operator with the Thomas--Fermi model.')
                g = this.CalculateSpeciesScreeningFactor( 1,Z,Ne,Z0,p,Lmax );
                return;
            elseif (useScreening == 5 && ~isfield(aTFDFT, speciesName)) ||...
                    (useScreening == 6 && ~isfield(FDFT, speciesName))
                %data missing from DFT or TFDFT, change to Boltz with Thomas--Fermi
                warning('no data on chosen species. Using Boltzmann operator with the Thomas--Fermi model.')
                g = this.CalculateSpeciesScreeningFactor( 4,Z,Ne,Z0,p,Lmax );
                return;
            end
            %__________________________________________________________________________
            %now calculate the screeningfactor as a function of Lvec = (1:Lmax)'
            switch useScreening
                case 0 %no screening
                    g = zeros(Lmax, length(p));
                case 1 %FP, TF
                    a_lengthscale = aTF_Krillov(Z,Ne);
                    g = ones(Lmax,1)*calcg_TFDFT_p(p,a_lengthscale);
                case 2 %FP, TFDFT
                    a_lengthscale = aTFDFT.(speciesName);
                    g = ones(Lmax,1)*calcg_TFDFT_p(p,a_lengthscale);
                case 3 %FP, DFT
                    % it is better to interpolate in log p space because we know that g
                    % increases logarithmically at large p. Otherwise the logarithm
                    % will cause g to decrease when extrapolating
                    g = ones(Lmax,1)*pchip(log(gDFT.p), gDFT.(speciesName), log(p));
                case 4 %Boltz, TF
                    %calculate the screening factor from scratch
                    % this takes a while (~8 s for 1000 momentum points) if want to use
                    % the same species many times consider saving and loading it?
                    a_lengthscale = aTF_Krillov(Z,Ne);
                    [gToInterp, pToInterp] = CalcgBoltzmann( @(p)calcF_TFDFT_p(p,a_lengthscale),Lmax, Z0, Ne );
                    g = pchip(log(pToInterp), gToInterp, log(p));
                case 5 %Boltz, TFDFT
                    %note if the matrix has dim [L, length(p)], the interpolation
                    %does what we want!
                    if isfield(gBoltzTFDFT, speciesName) && Lmax<=200
                        gAllL = pchip(log(gBoltzTFDFT.p), gBoltzTFDFT.(speciesName), log(p));
                        g = gAllL(1:Lmax, :);
                    else %calculate directly from F
                        a_lengthscale = aTFDFT.(speciesName);
                        [gToInterp, pToInterp] = CalcgBoltzmann( @(p)calcF_TFDFT_p(p,a_lengthscale),Lmax, Z0, Ne );
                        g = pchip(log(pToInterp), gToInterp, log(p));
                    end
                case 6 %Boltz, DFT
                    %note if the matrix has dim [L, length(p)], the interpolation
                    %does what we want!
                    if isfield(gBoltzDFT, speciesName)  && Lmax<=200
                        gAllL = pchip(log(gBoltzDFT.p), gBoltzDFT.(speciesName), log(p));
                        g = gAllL(1:Lmax, :);
                    else %calculate directly from F
                        calcF_DFT = pchip(log(FDFT.p), FDFT.(speciesName));
                        [gToInterp, pToInterp] = CalcgBoltzmann( @(p)ppval(calcF_DFT,log(p)), Lmax, Z0, Ne );
                        g = pchip(log(pToInterp), gToInterp, log(p));
                    end

                otherwise
                    warning('invalid useSpecies option chosen. Must be an integer 0-6.\n')
                    g = zeros(Lmax, length(p));
            end
        end

        function screeningFactor  = CalculateScreeningFactor(this, species, p, Lmax, LnLambda, useScreening, iteration)
            %CALCULATESCREENINGFACTOR calculates the enhancement of the electron-ion
            % deflection frequency from screening effects
            %   INPUTS:
            %     species has the fields
            %     nj, Z0NetCharge, ZAtomicNumber, times:
            %     nj(jspecies, times) - matrix containing number density of each
            %       species as a function of time (m^{-3})
            %     ZAtomicNumber - atomic number of each species. row vector, not a
            %     function of time
            %     Z0NetCharge - net charge for each species. row vector, not a
            %     function of time
            %     times - time steps. If constant, just set i to 1.
            %     - p = normalized momentum gamma*m*v/(m*c) = delta*y
            %     - LnLambda the coulomb logarithm. It is now a
            %       (Nxi-1, Ny ) matrix! (p-dependent)
            %     - useScreening: which model to use
            %       - 0: Ignore screening (complete screening)
            %       - 1: Fokker--Planck collision operator, TF model (works for all species)
            %       - 2: Fokker--Planck collision operator, TF-DFT model (works for some
            %			ions)
            %       - 3: Fokker--Planck collision operator, full DFT model (works for some
            %			ions)
            %       - 4: Boltzmann collision operator, TF model
            %       - 5: Boltzmann collision operator, TF-DFT model
            %       - 6: Boltzmann collision operator, full DFT model
            %       - 7: no screening = full penetration
            %           ADD WARNING IF P IS TOO BIG!!!
            %   OUTPUT: screeningFactor is a matrix of size (Nxi, Ny), that
            %   describes the enhancement of the deflection frequcency for the
            %   each Legendre mode and all y values on the state.

            nj_t = species.nj(:,iteration); %density of each species
            Z0NetCharge = species.Z0NetCharge;
            ZAtomicNumber = species.ZAtomicNumber;
            Ne = ZAtomicNumber-Z0NetCharge; %number of bound elecrons

            zWeight = (Z0NetCharge.^2)' * nj_t; %to divide in when calculating the screening factor
            screeningFactor = 1;

            if useScreening == 0 %full screening
                screeningFactor = ones(Lmax, length(p));
            elseif useScreening ==7 %complete penetration
                screeningFactor = 1/(zWeight) * ((ZAtomicNumber.^2)' * nj_t)...
                    *ones(Lmax, length(p));
            else %use our screening models
                for ispecies = 1:numel(nj_t)%for every species
                    % the contribution to screening for each species:
                    g = this.CalculateSpeciesScreeningFactor( useScreening,...
                        ZAtomicNumber(ispecies),Ne(ispecies),Z0NetCharge(ispecies),p,Lmax );

                    screeningFactor = screeningFactor + ...
                        1/(zWeight) * nj_t(ispecies) .* (g./LnLambda);
                end
            end
        end


       function [ boltzmanng,p ] = CalcgBoltzmann(~, F,Lmax, Z0, Ne )
            %CALCG_BOLTZMANN calculates the Legendre-mode-dependent correction to the
            %Fokker-Planck deflection frequency neglecting screening.
            %   input:
            %       - F(p): a function handle (if numeric, pchip(F,p) works
            %       - k =  p/mc (vector)
            %       - Z0, Ne: net charge and number of bound electrons
            %       - Lmax: maximum legendre mode

            % % The integrals are split at \theta_cutoff
            % because the legendres behave badly for too small arguments

            %definitions:
            cutoff = 1e-5; % before this P_L(1-2x^2)/x^2 = L(L-1) is a very
            % good approximation which is more stable than the polynomial evaluation
            Lvec = (1:Lmax)'; %vector of legendre mode numbers

            nMomentumPoints = 1e3;
            p = logspace(-4,3,nMomentumPoints);

            %the Z0 (no screening) integral (Z0^2*log(cutoff) below is part of this integral):
            BfactorIntegrand1 = @(x)Z0^2/x^3*...
                (1-LegendrePolynomials(Lmax,1-2*x^2))...
                *((1-x^2)*p.^2+1);
            BfactorNoScreening = integral(@(x)BfactorIntegrand1(x), cutoff, 1, 'arrayvalued', true);

            %The non-Z_0 parts, due to screening. Do not diverge, so we extend the
            % integral to 0, to make it independent of lnLambda.

            %small arguments:
            BfactorIntegrandScreening1 = @(x)...
                ( (2*Z0.*(Ne-F(x*p)) + (Ne-F(x*p)).^2)./(x*ones(size(p))) );
            %larger arguments:
            BfactorIntegrandScreening2 = @(x)1/x^2*(1-LegendrePolynomials(Lmax,1-2*x^2))*...
                ((2*Z0.*(Ne-F(x*p)) + (Ne-F(x*p)).^2)/x.*...
                ((1-x^2)*p.^2+1));
            BfactorScreeningPart1 = integral(@(x)BfactorIntegrandScreening1(x), 0, cutoff, 'arrayvalued', true);
            BfactorScreeningPart2 = integral(@(x)BfactorIntegrandScreening2(x), cutoff, 1, 'arrayvalued', true);

            % add it up
            boltzmanng = Z0^2*log(cutoff) + ...
                1./((Lvec.*(1+Lvec))*(p.^2+1)) .* BfactorNoScreening(2:end,:) +...
                ones(size(Lvec))* BfactorScreeningPart1 +...
                1./((Lvec.*(1+Lvec))*(p.^2+1)) .* BfactorScreeningPart2(2:end,:);%remove zeroth L mode where used LegendrePolynomials

        end

    end

end
