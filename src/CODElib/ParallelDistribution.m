function [f, y] = ParallelDistribution(dist)
% Calculates the parallel distribution (i.e. the xi=1 part). Wrapper for
% SumLegModesAtXi1.
[posF, negF] = SumLegModesAtXi1(dist.f, dist.momentumGrid.Ny, dist.momentumGrid.Nxi);

y = [-flip(dist.momentumGrid.y) dist.momentumGrid.y];
f = [flip(negF); posF];
end