classdef RosenbluthSource < Source
    %ROSENBLUTH Rosenbluth Putvinski like source.
    % Can either use number of runaways or number of fast particles
    
    methods
        function this = RosenbluthSource(state,eqSettings)
                this@Source(state,eqSettings)
        end
        
        function sourceVector = getSourceVec(this,f,iteration)
            
            deltaRef        = this.state.reference.deltaRef;
            lnLambdaRef     = this.state.reference.lnLambdaRef;
            
            EOverEc         = this.state.physicalParams.EOverEc(iteration);
            delta           = this.state.physicalParams.deltas(iteration);
            nBar            = this.state.physicalParams.nBars(iteration);
            neTotalOverneFree = this.state.physicalParams.neTotalOverneFree(iteration);
            veBar3          = this.state.physicalParams.veBars3(iteration);
            
            y               = this.state.momentumGrid.y;
            yWeights        = this.state.momentumGrid.yWeights;
            x               = this.state.momentumGrid.x;
            gamma           = this.state.momentumGrid.gamma;
            Ny              = this.state.reference.momentumGrid.Ny;
            Nxi             = this.state.reference.momentumGrid.Nxi;
            matrixSize      = this.state.momentumGrid.matrixSize;
                       
            sourceMode      = this.eqSettings.sourceMode;
            
            w               = deltaRef*y;
            
            if abs(EOverEc) > 1
                gammaC = sqrt(abs(EOverEc)/(abs(EOverEc) - 1));
            else
                gammaC = inf;
            end
            wMinForSource   = sqrt(gammaC^2-1);
            nr_mask = w > wMinForSource;
            %before main loop
            if sourceMode == 2
                [~,tail_mask] = DetermineFastParticleRegion(this,x,y,f,...
                                                            nBar,veBar3,delta,deltaRef,nr_mask); %switched 0 to nr_mask
            end


            mask0 = deltaRef*y>delta*this.eqSettings.yCutSource;
            if this.eqSettings.sourceMode == 2
                if this.eqSettings.fastParticleDefinition == 4 && iteration == 1  
                    mask = mask0 | tail_mask;
                else
                    mask = tail_mask;
                end
            else
                mask = mask0;
            end
            
            %Preliminaries for matrix
            RosenbluthSourceVector0 = zeros(matrixSize,1);
            if this.eqSettings.sourceMode == 1 || this.eqSettings.sourceMode == 2
                xi2 = -sign(EOverEc) * w./(1+gamma);    
                prefactor = deltaRef^5 * 3*pi/(lnLambdaRef*16);
                RPSourceFunc = prefactor./(y.*gamma.*(gamma-1).*(gamma-1));
                RPSourceFunc(1) = 0;
                if this.state.momentumGrid.yMaxBoundaryCondition ~= 3
                    RPSourceFunc(end) = 0;
                end
                PLMat = LegendrePolynomials(Nxi-1,xi2);
                for L=0:(Nxi-1)
                    indices = L*Ny+(1:Ny);
                    RosenbluthSourceVector0(indices) = (2*L+1)/2*PLMat(L+1,:) .* RPSourceFunc;
                end
            end
            
            %in main loop
            
            RosenbluthSourceVector = nBar * neTotalOverneFree * RosenbluthSourceVector0 .* repmat(mask(:),Nxi,1);                    
            f0 = f(1:Ny);   
            
            pc_RP        = deltaRef * y(find(RosenbluthSourceVector,1));
            if ~isempty(pc_RP)
                    gammaC_RP    = sqrt(1+pc_RP^2);
                    pc_RP_star   = sqrt(4*gammaC_RP*(gammaC_RP-1));
                    nr_RP        = 4/sqrt(pi) * (yWeights.*y.^2.*(deltaRef*y > pc_RP_star))*f0;
                    sourceVector = nr_RP*RosenbluthSourceVector;
            else
                    sourceVector = zeros(matrixSize,1);
            end
        end


        function [idTail,tail_mask] = DetermineFastParticleRegion(this,x,y,f,...
                                                    nBar,veBar3,delta,deltaRef,nr_mask)
            switch this.eqSettings.fastParticleDefinition
                case 0 %Do not calculate fast particles
                    idTail = [];
                case 1 %Based on the beginning of the tail 
                    [fAtXi1,~] = SumLegModesAtXi1(f,Ny,Nxi);
                    maxw = nBar/veBar3*exp(-y.^2/veBar2);
                    quota = fAtXi1'./maxw;
                    idTail = find(quota>this.eqSettings.tailThreshold,1);
                case 2 %Based on relative speed 
                    idTail = find(deltaRef*x > delta*this.eqSettings.relativeSpeedThreshold,1);
                case 3 %Based on absolute speed
                    idTail = find(deltaRef*x > this.eqSettings.absoluteSpeedThreshold,1);
                case 4 %Combination                    
                    idRE = find(nr_mask,1);
                    idRel = find(deltaRef*x > delta*this.eqSettings.relativeSpeedThreshold,1);
                    idAbs = find(deltaRef*x > this.eqSettings.absoluteSpeedThreshold,1);
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
    end
end
