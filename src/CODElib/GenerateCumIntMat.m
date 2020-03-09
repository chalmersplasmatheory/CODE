function [m,xiGrid] = GenerateCumIntMat(Nxi,nPoints)
    %Cumulatively integrate legendre modes over a xi grid

    xiGrid = linspace(0,1,nPoints);
    dXi = xiGrid(2)-xiGrid(1);
    legPols = LegendrePolynomials(Nxi-1,xiGrid);
    m = dXi*fliplr(cumtrapz(fliplr(legPols),2));        
end