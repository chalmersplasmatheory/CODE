function n = CalculateDensity(f,Ny,Nxi,y,yWeights,nRef)
%CALCULATEDENSITY calculates density n as integral over momentum space. The
%density is output in SI.

% %%%%%%%%%%%%%%%%
% Input:
% %%%%%%%%%%%%%%%%
% f - distribution function normalized as in CODE
% y - momentum grid
% yWeights - weights for integrating y
% nRef - reference density
% Nxi - Legendre Modes
% Ny - points in y;

% Written by Albert Johansson 2019
    y = y(:)';
    yWeights = yWeights(:)';
    f = reshape(f,Ny,Nxi);

    preFac = 2 * nRef / sqrt(pi);

    yParts = (yWeights.*y.^2) * (f);

    XiIntegrand = @(x) LegendrePolynomials(Nxi-1,x)';
    XiParts = integral(XiIntegrand,-1,1,'ArrayValued',true);

    n = preFac * yParts * XiParts';
end

