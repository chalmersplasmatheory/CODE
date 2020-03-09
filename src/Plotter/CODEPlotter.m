classdef CODEPlotter
    %CODEPLOTTER Static clacontaining functions for plotting results
    %from CODE runs or plotting mid run.
    
    properties
        figureOffset =1 
    end
    
    methods (Static)
        function  figureHandle = contourDistribution(distr, figureOffset)
            %contourDistribution plots the distribution on a countour graph
			%internal figure number 1.
            % distr - object of Distribution to be plotted
            Nxi = distr.momentumGrid.Nxi;
            w = distr.momentumGrid.y*distr.momentumGrid.reference.deltaRef;
            EOverEc = distr.EOverEc;
            if abs(EOverEc) < 1
                wc = Inf;
            else
                wc = 1/sqrt(abs(distr.EOverEc)-1);
            end
            if nargin < 2
                figureOffset = 0;
            end
			figureHandle = figure(1+figureOffset);
			clf
			numContours = 15;
			% Maximum momenutm which can be fitted in a square from circular coordinates with maximum radius
			% in contour plot without having to extrapolate, rounding sqrt(2) up slightly:
            % also remove w where distribution is not of any importance
            % size
			Nw = max(size(w));
			wMax = max(w(distr.f(1:Nw)<1e-15))/1.5;
            % Number of grid points to use in making the contour plots:
			NPlot = 201 ;
			% It's best to choose NPlot odd to avoid the singularity at the origin.
			
			wPar1D = linspace(-wMax, 0.5*wMax, NPlot);
			wPerp1D =linspace(-wMax, wMax, NPlot);
			[wPar2D, wPerp2D] = meshgrid(wPar1D, wPerp1D);
			w2D = sqrt(wPar2D.*wPar2D + wPerp2D.*wPerp2D);
			xi2D = wPar2D ./ w2D;
			
			w2D = reshape(w2D,[NPlot*NPlot,1]);
			xi2D = reshape(xi2D,[NPlot*NPlot,1]);
			
			f = zeros(NPlot*NPlot,1);
			for L=0:(Nxi-1)
				fSlice = distr.f(L*Nw+(1:Nw));
				
				% Use Matlab's built-in function to evaluate the Legendre
				% polynomials on the points we need for the contour plot:
				LegendreMatrix = legendre(L, xi2D);
				% Matlab's 'legendre' function also returns a bunch of other
				% uninteresting numbers besides the Legendre polynomials, so 
				% keep only the values of the Legendre polynomials:
				Legendres = LegendreMatrix(1,:)';
				
				% Interpolate from the regular y grid to the points we need in
				% the Cartesian plane for the contour plot:
				f = f + Legendres .* interp1(w, fSlice, w2D);
			end
			
			f = reshape(f, [NPlot, NPlot]);
			
			toPlot = log10(abs(f));
			minContour = max([-15, min(min(toPlot))]);
			contours = linspace(minContour, 0, numContours);			
			contourf(-wPar2D, wPerp2D, toPlot, contours)
			colorbar
			xlabel('w_{||}')
			ylabel('w_{\perp}')
			title('log_{10}|F|')
			%%%%%%%%%%%%%%%%%%%% plot different runaway regions %%%%%%%%%%%%%%%%%%%%%
			hold on;
			%Isotropic
			th = linspace(0,2*pi,1000);
			circPar = cos(th);
			circPerp = sin(th);
			h1 = plot(wc*circPar,wc*circPerp,'--r','LineWidth',1);
			
			%Region in Smith et al. (valid non-realtivistically) 
			xi = circPar;
			wSep = sqrt(1./abs(EOverEc))./sqrt((xi+1)/2);
			wSPar = xi.*wSep;
			wSPerp = wSep.*sqrt(1-xi.*xi);
			wSPerp(wSPar<-wMax) = [];
			wSPar(wSPar<-wMax) = [];
			wSPar (wSPerp>max(wPerp2D(:))) = [];
			wSPerp(wSPerp>max(wPerp2D(:))) = [];
			plot(wSPar,wSPerp,'-.c','LineWidth',2);
			h3 = plot(wSPar,-wSPerp,'-.c','LineWidth',2);
			
			%Xi-dependent
			xi = circPar(circPar>0);
			arg = xi*EOverEc-1;
			xi(arg<=0) = [];
			arg(arg<=0) = [];
			wcX = 1./(sqrt(arg));        
			wcPar = xi.*wcX;
			wcPerp = wcX.*sqrt(1-xi.*xi);
			wcPar(wcPerp>max(wPerp2D(:))) = [];
			wcPerp(wcPerp>max(wPerp2D(:))) = [];
			plot(wcPar,wcPerp,'--y','LineWidth',2);
			h2 = plot(wcPar,-wcPerp,'--y','LineWidth',2);
			
			
			
			legend([h1,h2,h3],'w_c (isotropic)','w_c(\xi)','x_{sep} (Smith)')
			
			axis equal 
            drawnow;
        end
        
        function Plotp2f(hDistEvo,distr, FMinForPlot, yMaxForPlot)
            %hDistEvo - plot handle
            % distr - object of distribution class
            % FMinForPlot - minimum distr value for plot, optional
            % yMaxForPlot - max momentum for plot, optional
            y = distr.momentumGrid.y;
            deltaRef = distr.momentumGrid.reference.deltaRef;
%             EOverEc = distr.EOverEc;
            fMinus1 = distr.f;
            if nargin <3
                FMinForPlot = 1e-18;
            end
            if nargin < 4
               yMaxForPlot = 50; 
            end
            Ny = numel(y);
            cla(hDistEvo,'reset'); 
            semilogy(deltaRef*y,(deltaRef*y').^2.*fMinus1(1:Ny),'--b','linewidth',2)
            hold on
            FDistPlot = (deltaRef*y').^2.*fMinus1(1:Ny);
            hDistLine = semilogy(deltaRef*y,FDistPlot,'k','linewidth',3);
            hDistPcLine = semilogy([0 0],FMinForPlot*[1 1],':k','linewidth',2);
            xlim(deltaRef*[0,yMaxForPlot]) 

            hP_thresh = min(y(end)*deltaRef/2,0.1);
            hI_thresh = find(deltaRef*y>hP_thresh,1);
            %ylim([FMinForPlot,max(1e10*FMinForPlot, 2*max((deltaRef*y(hI_thresh:Ny)').^2.*fMinus1(hI_thresh:Ny)))])
            ylim([FMinForPlot, 20])
            
            set(gca,'fontsize',16,'ticklabelinterpreter','latex','linewidth',2)
            xlabel('$p/m_e c$','fontsize',18,'interpreter','latex')
            ylabel('$p^2 f_0$','fontsize',18,'interpreter','latex')
            leg = legend('Initial value','Distribution','$p_c$ (Connor-Hastie)');
            set(leg,'interpreter','latex','box','off','fontsize',18) 
        end
        
        
        function PlotInstantaneousDist(plotHandle,distr,FMinForPlot, yMaxForPlot)
            clf(plotHandle); 
            hDistEvo = axes;
            y = distr.momentumGrid.y;
            soln = distr.f;
            Ny = distr.momentumGrid.Ny;
            Nxi = distr.momentumGrid.Nxi;
            yCrit = 1/sqrt(abs(distr.EOverEc)-1)/distr.momentumGrid.reference.deltaRef;
            nBar = distr.nBar;
            veBar3 = distr.veBar3;
            time = distr.time;

            if nargin <3
                FMinForPlot = 1e-18;
            end
            if nargin < 4
               yMaxForPlot = 50; 
            end
            
            colors = {'b','r','g','c'};
            colorsP = {[1,0.8,0],'m','b','r'};
            warning off MATLAB:Axes:NegativeDataInLogAxis

            %Positive and negative here has been reversed since the code was
            %written. Negative X here thus means the direction of the runaway
            %tail.
            cla(hDistEvo);

            [FNegativeX,FPositiveX] = SumLegModesAtXi1(soln,Ny,Nxi);
            xBig = -[fliplr(-y), y]';
            F = [flipud(FNegativeX); FPositiveX];
            semilogy(hDistEvo,xBig, F, 'Color',colors{1});
            hold(hDistEvo,'on');        
            plotHandle = semilogy(hDistEvo,xBig, -F, ':','Color',colors{1},'LineWidth',2);
            set(get(get(plotHandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend        
%             if size(varargin) == 1 %If there is a primary-only-distribtuion
%                 [FNegativeX,FPositiveX] = SumLegModesAtXi1(varargin{1},Ny,Nxi);            
%                 FPrim = [flipud(FNegativeX); FPositiveX];
%                 semilogy(hDistEvo,xBig, FPrim,'--','Color',colorsP{1});                  
%                 plotHandle = semilogy(hDistEvo,xBig, -FPrim, ':','Color',colorsP{1},'LineWidth',2);
%                 set(get(get(plotHandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend        
%             end

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
            legend(hDistEvo,legT);
            title(hDistEvo,['Total distribution for v_{perp} = 0 at t = ',num2str(time),' reference collision times']);        
            xlabel(hDistEvo,'\gamma v_{||} / v_{th,ref}')
            ylabel(hDistEvo,'F')
    %         pause(0.2)
            drawnow;

        end

    end
end

