classdef ImprovedChiuHarveySource  < ExplicitOperator
    %ChiuHarvey45Source

    properties
        NyInterp = 1000;
    end
    
    %used in knock on matrix
    properties
        deltaRefUsed = 0
        yUsed = 0
        deltaUsed = 0
        EOverEcUsed = 0
        pCutSourceUsed = 0
        NxiUsed = 0
        sourceModeUsed = 0
    end
    %used outside knock on matrix
    properties
        lnLambdaRefUsed = 0
        nBarUsed = 0
        neTotalOverneFreeUsed = 0
    end
    
    methods
        function this = ImprovedChiuHarveySource(state,eqSettings)
            this@ExplicitOperator(state, eqSettings)
        end
    end

    methods

        function generateOperatorMatrix(this,iteration)
            
            lnLambdaRef         = this.state.reference.lnLambdaRef;
            min_delta           = min(this.state.physicalParams.deltas);
            delta               = this.state.physicalParams.deltas(iteration);
            pCutSource          = this.eqSettings.yCutSource*min_delta;
            %Knock-on source
            if this.eqSettings.sourceMode == 5 || this.eqSettings.sourceMode == 4 && (pCutSource < deltaRef*this.state.momentumGrid.yMax)

                % preFac = 0.25*3*sqrt(pi) * delta^3/lnLambda;
                % Ola comment 2015-02-16: I just divided this by an additional
                % factor of two, which appeared to be missing from the
                % expression. At some point double-check normalization, perhaps.
                % At least now source 4 and 5 produce nearly the same growth
                % rate, rather than being off by a factor 2.
                preFac = (3*sqrt(pi)/8) * delta^3/lnLambdaRef;
                
                if this.deltaRefUsed == this.state.reference.deltaRef && ...
                            isequal(this.yUsed,this.state.momentumGrid.y) && ...
                            this.deltaUsed == this.state.physicalParams.deltas(iteration) && ...
                            this.EOverEcUsed == this.state.physicalParams.EOverEc(iteration) && ...
                            this.pCutSourceUsed == min_delta*this.eqSettings.yCutSource && ...
                            this.NxiUsed == this.state.momentumGrid.Nxi && ...
                            this.sourceModeUsed == this.eqSettings.sourceMode
                    koMat = this.operatorMatrix...
                        /this.deltaUsed^3*this.lnLambdaRefUsed/(3*sqrt(pi)/8)*preFac...
                        /this.nBarUsed/this.neTotalOverneFreeUsed;
                else
                    koMat = preFac*this.GenerateKnockOnMatrix(pCutSource);
                end
            else
                sz = this.state.momentumGrid.Nxi*this.state.momentumGrid.Ny;
                koMat = sparse(sz,sz);
            end
            
            nBar = this.state.physicalParams.nBars(iteration);
            neTotalOverneFree = this.state.physicalParams.neTotalOverneFree(iteration);
            this.operatorMatrix = nBar*neTotalOverneFree*koMat;
            
            this.deltaRefUsed = this.state.reference.deltaRef;
            this.yUsed = this.state.momentumGrid.y;
            this.deltaUsed = delta;
            this.EOverEcUsed = this.state.physicalParams.EOverEc(iteration);
            this.pCutSourceUsed = this.eqSettings.yCutSource*min_delta;
            this.NxiUsed = this.state.momentumGrid.Nxi;
            this.sourceModeUsed = this.eqSettings.sourceMode;
            this.lnLambdaRefUsed = lnLambdaRef;
            this.nBarUsed = nBar;
            this.neTotalOverneFreeUsed = neTotalOverneFree;
        end

        function [M,M_source,M_sink1,M_sink2] = GenerateKnockOnMatrix(this, pm)
            disp('Knockon matrix builds')
            p                   = this.state.reference.deltaRef*this.state.momentumGrid.y;
            yWeights            = this.state.momentumGrid.yWeights;
            deltaRef            = this.state.reference.deltaRef;
            totalTime           = tic;

            numInterpPoints     = this.eqSettings.NyInterp;
            Ny                  = this.state.momentumGrid.Ny;
            Nxi                 = this.state.momentumGrid.Nxi;
            p                   = p(:);

            % overwrites the above loop stuff
            pm_Ind = find(p>pm,1);

            pWeights = deltaRef*yWeights(:);   

            pWeights(end) = pWeights(end-1);

            pMax     = p(end);
            gamma    = sqrt(1+p.^2);
            gamma_m  = gamma(pm_Ind);    

            P        = repmat(p,1,numInterpPoints);
            Gamma    = sqrt(1+P.^2);
            P_interp = pMax*ones(size(P));
            WS       = zeros(size(P));


            useReducedKnockOnSource = (this.eqSettings.sourceMode == 4); % Generalized Chiu-Harvey operator (neglects test-particle piece) on matrix form

            if useReducedKnockOnSource
                p_int_min = sqrt( (2*gamma-1).^2-1 );    
            else
                p_int_min = sqrt( (gamma+gamma(pm_Ind)-1).^2-1 );    
            end
            p_int_max = pMax*ones(size(p));


        %     calculateOnlyTestParticle = 0;
        %     if calculateOnlyTestParticle 
        %         p_int_min = sqrt( (gamma+gamma(pm_Ind)-1).^2-1 );    
        %         p_int_max = min(sqrt( (2*gamma-1).^2-1 ),pMax);  
        %         P_interp  = repmat(p_int_max,1,numInterpPoints); %zeros(size(P));
        %     end
            %gamma_int_max = gamma(end) - (gamma_m-1);
            %p_int_max = sqrt(gamma_int_max^2-1)*ones(size(p));

            I_lower      = pm_Ind;
            I_upper      = find(gamma>=(gamma(end) + 1 - gamma_m),1);
            parfor K = I_lower:I_upper         % Interpolating to a Gauss-Legendre
                                                % quadrature grid (using lgwt).
                if p_int_max(K) > p_int_min(K)
        %             if useReducedKnockOnSource
        %                 I_min = find(p>p_int_min(K),1);
        %                 ps        = linspace(p(I_min),p_int_max(K),numInterpPoints);
        %                 ws        = ps(2)-ps(1);
        %             else
                        [ps,ws]   = lgwt(numInterpPoints,p_int_min(K),p_int_max(K));
                        ps = fliplr(ps);
                        ws = fliplr(ws);
        %             end
                    P_interp(K,:) = ps;
                    WS(K,:)       = ws;
                end
            end
            Gamma_interp = sqrt(1 + P_interp.^2);

            costheta_s   = sqrt( (Gamma_interp+1).*(Gamma-1)./( (Gamma_interp-1).*(Gamma+1) ) );
            PL           = LegendrePolynomials(Nxi-1,costheta_s(:)');
            PL           = reshape(PL,Nxi,Ny,numInterpPoints);
            Sigma        = CalculateMollerCS(Gamma,Gamma_interp);

            ML0          = zeros(Ny,numInterpPoints);

            indices         = I_lower:I_upper;
            ML0(indices,:)  = WS(indices,:).*P_interp(indices,:).^3 ./(Gamma_interp(indices,:) ...
                               .*P(indices,:).*Gamma(indices,:)).* Sigma(indices,:);%...
        %   .*(Gamma(indices,:)>=gamma_m).*(Gamma(indices,:) <= (gamma(end) + 1 - gamma_m));
                                                % The common momentum 
                                                % dependence for all L.

            ML0(1,:) = 0;                       % Eliminating division by 
                                                % zero with the factor 1/P

        %     progressPhrases = {'Sorting variables alphabetically...', ...
        %                        'Converting numbers to hexadecimal...',...
        %                        'Counting switch statements in CODE...',...
        %                        'Trying a Lagrangian approach...',...
        %                        'Starting from scratch...',...
        %                        'Estimating chance of success...',...
        %                        'Deriving analytic estimate...'};


            MemorySaver = 2;                    % REVISE THIS ADVICE:
                                                % Use MemorySaver for ~half the memory 
                                                % consumption, ~10 times computation time
            if MemorySaver==0
        %         t_memSave = tic;
                I   = zeros((Ny-pm_Ind)*Ny*Nxi,1);  % The full I,J will represent Nxi solid
                J   = I;                            % blocks along the diagonal.
                for K = pm_Ind:Ny  
                    % This is indeed a pretty funky way of indexing the matrix. 
                    % I am building it row-by-row for all Legendre modes simultaneously.
                    indices      = (K-pm_Ind)*Ny*Nxi + (1:(Ny*Nxi));
                    Lmat         = repmat(0:(Nxi-1),Ny,1);
                    Imat         = K*ones(Ny,Nxi) + Ny*Lmat;
                    Jmat         = repmat((1:Ny)',1,Nxi) + Ny*Lmat;
                    I(indices)   = Imat(:);
                    J(indices)   = Jmat(:);
                end
                S = cell(Ny,1);
                parfor K = I_lower:Ny
                    MLt = zeros(Nxi,Ny);
                    if p_int_max(K) > p_int_min(K)
                        interpMatrix    = this.BuildInterpolationMatrix(p,P_interp(K,:));
                        ML0t            = repmat(squeeze(ML0(K,:)),Nxi,1);
                        MLt(:,:)        = (ML0t.*squeeze(PL(:,K,:))) * interpMatrix;
                        MLt = MLt';
                    end
                    S{K} = MLt(:);    
                end
                S = cell2mat(S);
                clear('PL')
                pause(1e-6)
                M_source = sparse(I,J,S,Ny*Nxi,Ny*Nxi); 
                clear('S','I','J')
        %         toc(t_memSave);
            elseif MemorySaver==1
        %         t_memSave = tic;
                M_source = spalloc(Ny*Nxi,Ny*Nxi,Ny^2*Nxi/2);
                Lmat     = Ny*repmat(0:(Nxi-1),Ny,1);
                oneMat   = ones(Ny,Nxi);
                NMat     = repmat((1:Ny)',1,Nxi);
                Jmat     = NMat + Lmat;
                J        = Jmat(:);
                counter = 0;
                fprintf('Progress on generating the knock-on matrix:\n') 
                for K = 1:Ny
                    Imat    = K*oneMat + Lmat;
                    I       = Imat(:);
                    if p_int_max(K) > p_int_min(K)
                        interpMatrix    = this.BuildInterpolationMatrix(p,P_interp(K,:));
                        ML0t            = repmat(squeeze(ML0(K,:)),Nxi,1);
                        MLt             = (ML0t.*squeeze(PL(:,K,:))) * interpMatrix;
                        MLt             = MLt';
                        M_source        = M_source + sparse(I,J,MLt(:),Ny*Nxi,Ny*Nxi);
                    end
                    if mod(K,round(Ny/10))==0
                        counter = counter+1;
                        if counter == 6
                            fprintf('\n')
                        end
                        if counter == 10
                            fprintf('100%%. Done!\n')
                        else
                            fprintf('%0.0f%%... ',10*counter)
                        end
                    end
                end
                clear('PL')
        %         toc(t_memSave)
            elseif MemorySaver==2 % Works on up to blockSize number of threads in parallel, 
                                  % successively builds the matrix for blockSize
                                  % number of K values at a time.
        %         t_memSave = tic;
                blockSize = 64;

                Lmat     = Ny*repmat(0:(Nxi-1),Ny,1);
                oneMat   = ones(Ny,Nxi);
                NMat     = repmat((1:Ny)',1,Nxi);
                Jmat     = NMat + Lmat;
                J0       = Jmat(:);
                M_source = spalloc(Ny*Nxi,Ny*Nxi,round((I_upper-I_lower)^2*Nxi/2));

                maxBlockNumber = floor(Ny/blockSize);
                remainderSize = mod(Ny,blockSize);
                for blockNumber = 1:maxBlockNumber
                    I   = cell(blockSize,1);
                    S   = cell(blockSize,1);

                    P_interp_array  = cell(blockSize,1);
                    ML0_array       = cell(blockSize,1);
                    PL_array        = cell(blockSize,1);
                    p_int_maxes     = p_int_max((1:blockSize) + blockSize*(blockNumber-1));
                    p_int_mins      = p_int_min((1:blockSize) + blockSize*(blockNumber-1));
                    for R = 1:blockSize;
                        K = R + blockSize*(blockNumber-1);
                        P_interp_array{R} = P_interp(K,:);
                        ML0_array{R}      = ML0(K,:);
                        PL_array{R}       = squeeze(PL(:,K,:));
                    end

                    parfor R = 1:blockSize
        %                 sprintf('R=%d, ',R)
                        K = R + blockSize*(blockNumber-1);
                        Imat    = K*oneMat + Lmat;
                        I{R}    = Imat(:);
                        S{R}    = zeros(size(I{R}));
                        if p_int_maxes(R) > p_int_mins(R)
                            interpMatrix    = this.BuildInterpolationMatrix(p,P_interp_array{R});
                            ML0t            = repmat(squeeze(ML0_array{R}),Nxi,1);
                            MLt             = ((ML0t.*PL_array{R}) * interpMatrix)';
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
                    parfor R = 1:remainderSize
                        K = R + blockSize*maxBlockNumber;
                        Imat    = K*oneMat + Lmat;
                        I{R}    = Imat(:);
                        S{R}    = zeros(size(I{R}));
                        if p_int_max(K) > p_int_min(K)
                            interpMatrix    = this.BuildInterpolationMatrix(p,P_interp(K,:));
                            ML0t            = repmat(squeeze(ML0(K,:)),Nxi,1);
                            MLt             = ((ML0t.*squeeze(PL(:,K,:))) * interpMatrix)';
                            S{R} = MLt(:);
                        end
                    end
                    I = cell2mat(I);
                    S = cell2mat(S);
                    J = repmat(J0,remainderSize,1);
                    M_source = M_source + sparse(I,J,S,Ny*Nxi,Ny*Nxi);
                end
        %         toc(t_memSave)
                clear('S,I,J,MLt,PL,ML0t')
            end


            pause(1e-6)   


        %     sinkFunc1 = -0.5*((p.^2.*pWeights)'*M_source(1:Ny,1:Ny))./(p.^2.*pWeights)';
        %     sinkFunc1(isnan(sinkFunc1)) = 0;
        %     sinkFunc1(isinf(sinkFunc1)) = 0;
        %     sinkFunc2 = -((p.^3.*pWeights)'*M_source(Ny+1:2*Ny,Ny+1:2*Ny))./(p.^3.*pWeights)';
        %     sinkFunc2(isnan(sinkFunc2)) = 0;
        %     sinkFunc2(isinf(sinkFunc2)) = 0;


            c = 0; %shifting the cut-off energy for the sink term by a grid point
            gamma_m  = c*gamma(pm_Ind-1) + (1-c)*gamma(pm_Ind);



            sigma0   = CalculateSigmaIntegral(gamma,gamma_m);


            if useReducedKnockOnSource
                M_sink1 = 0*M_source;
            else
                M0_sink1Func        = -0.5*p./gamma .* (gamma>=(2*gamma_m-1)).*sigma0;
                S_sink1             = repmat(M0_sink1Func(:),Nxi,1);

        %     S_sink1(1:Ny)       = sinkFunc1;
        %     S_sink1(Ny+1:2*Ny)  = sinkFunc2;

                M_sink1             = spdiags(S_sink1,0,Ny*Nxi,Ny*Nxi);
            end


            a               = 2/deltaRef;          % scaling the bulk particle-sink, 
                                                % 1/delta corresponds to being 
                                                % proportional to the background Maxwellian.                  
            deltaSinkFunc   = 4*a^3/sqrt(pi) * (5/2 - (a*p).^2).*exp(-(a*p).^2); 
            normFactor2     = (p.^2.*pWeights)'*deltaSinkFunc;
            deltaSinkFunc   = deltaSinkFunc / normFactor2;
                                                % this function has zero (nonrelativistic)
                                                % energy moment and unity density moment.



            M0_sink2        = -0.5*deltaSinkFunc ...
                * ( (gamma' >= (2*gamma_m-1)) .* p'.^3 ./ gamma' .* sigma0' .* pWeights');


        %     sinkFunc2 = -0.5*((p.^2.*pWeights)'*M_source(1:Ny,1:Ny));
        %     M0_sink2 = deltaSinkFunc * sinkFunc2(:)';


            [J0,I0] = meshgrid(1:Ny);           
            I0      = I0(:);                    
            J0      = J0(:);
            S_sink2 = M0_sink2(:);





             M_sink2     = sparse([I0;(Ny*Nxi);round(Ny*Nxi/2)],[J0;round(Ny*Nxi/2);(Ny*Nxi)],[S_sink2;0;0]);
        %    M_sink2     = sparse(I0,J0,S_sink2,Ny*Nxi+1,Ny*Nxi+1);%[I0;(Ny*Nxi);round(Ny*Nxi/2)],[J0;round(Ny*Nxi/2);(Ny*Nxi)],[S_sink2;0;0]);

                                                    % Adding a zero-element on the last row
                                                    % and column to get the correct matrix
                                                    % size. Any less ugly way to do it??
            M = M_source + M_sink1 + M_sink2;


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



            fprintf('Calculated knock-on matrix in %0.1f seconds.\n',toc(totalTime))
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
        end

        function K = BuildInterpolationMatrix(~,x,xp)
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

