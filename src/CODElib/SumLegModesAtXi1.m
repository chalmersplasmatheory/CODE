function [posF,negF] = SumLegModesAtXi1(f,Ny,Nxi)
    %Calculates the total distribution at xi=1 (v_perp=0) by summing
    %over all the Legendre modes of the distribution

    % %%%%%%%%%%%%%%
    % Input:
    % %%%%%%%%%%%%%%
    % f - distribution function with normalization as in CODE
    % and Nxi legandre modes projections. first Ny points are for first
    % legandre mode.
    % Nxi - Number of Legandre Modes
    % Ny - Number of radial grid points

    distMatrix = reshape(f,Ny,Nxi);
    plusMin = ones(1,Nxi);
    plusMin(2:2:end) = -1;
    posF = sum(distMatrix*diag(plusMin),2); %Sum over all the legendre modes, with alternating sign. Creates a vector of values at the grid points
    negF = sum(distMatrix,2);
end
