function flowRate = GetFlowThroughGridBoundary(dist, yBoundary)
	% Calculates flow through the boundary y_b. Does not take into
	% account for several effects, such as screening. Should not affect
	% results much due to the electric-field term dominating for large
    % yBoundary.
    
    Ny = dist.momentumGrid.Ny;
    ddy = dist.momentumGrid.ddy;
    
    y = dist.momentumGrid.y;
    y2 = dist.momentumGrid.y2;
    x = dist.momentumGrid.x;
    x2 = dist.momentumGrid.x2;
    
    deltaRef = dist.deltaRef;
    nBar = dist.nBars;
    
    Fn = @(n)  dist.f(Ny*n + 1:Ny*(n+1));
	F0 = Fn(0);
	F1 = Fn(1);
	F2 = Fn(2);

	EHat = dist.EHats;
	BHat = dist.BHatRef;
	nueeBar = dist.nueeBars;
 	veBar = dist.veBars;

	dF0dy = (ddy*F0)';
	yb = find(y >= yBoundary, 1);
	veBar2 = veBar*veBar;
	expx2 = exp(-x2/veBar2);
	%         erfs = erf(x/veBar);
    erfs = erf(x/veBar);
	psi = (erfs-2/(veBar*sqrt(pi))*x.*expx2) * veBar2 ./ (2*x2);

    includeSynchrotronLoss = dist.B ~= 0;
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