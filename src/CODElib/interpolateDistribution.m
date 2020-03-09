function f = interpolateDistribution(dist,MG)
% Interpolates each legandre mode function from distribution for itself with the momentumGrid in query.
% The distribution in Distribution is interpreted to be zero above its maximum momentum, and in legandre modes higher than its Nxi.
% Input:
% dist - Distribution object, from which the distribution should be interpolated.
% MG - MomentumGrid object, to which grid the output distribution is interpolated.

if MG.y(end) > dist.momentumGrid.y(end)
    f = zeros(MG.matrixSize,1);
    for i = 0:min(MG.Nxi,dist.momentumGrid.Nxi)-1
        f(i*MG.Ny+1:(i+1)*MG.Ny) = interp1([dist.momentumGrid.y, MG.y(end)],[dist.f(i*dist.momentumGrid.Ny+1:(i+1)*dist.momentumGrid.Ny), 0],MG.y);
    end
else
    f = zeros(MG.matrixSize,1);
    for i = 0:min(MG.Nxi,dist.momentumGrid.Nxi)-1
        f(i*MG.Ny+1:(i+1)*MG.Ny) = interp1(dist.momentumGrid.y,dist.f(i*dist.momentumGrid.Ny+1:(i+1)*dist.momentumGrid.Ny),MG.y);
    end
end
end
