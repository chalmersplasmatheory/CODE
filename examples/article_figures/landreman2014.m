% This is an example script to recreate some of the figures from the
% first CODE article.
% I would recommend reading it to better understand this script.

%% Figure 1

% This figure uses Ehat and delta as units. For this to work with the
% latest version of CODE, we will have to convert them into SI units.

n = 1e19; % [m^-3]
E = 0.065; % [V/m]
T = 2.5550e+03; % [eV]
Z = 1;

% Since we will not change T och n during the simulation, we set the
% reference parameters to the same thing. The reference parameters
% determines how the grid is normalized.
S = State(T, n);
S.physicalParams.setPhysicalParameters(T, n, Z, E);

% The unit conversions can be confirmed by checking S.physicalParams.EHats,
% and S.physicalParams.delta.

% The same resolution paramters are used as in the article.
yMax = 20;
Ny = 300;
Nxi = 20;
S.momentumGrid.setResolution('yMax', yMax, 'Ny', Ny, 'Nxi', Nxi);

% CODE is very efficient when the physical paramters aren't changed during
% the simulation. Which is why we can set tMax very high here.
tMax = 500;
S.timeGrid.setResolution('tMax',tMax,'dt',1);

es = EquationSettings(S);

% There are two solvers in CODE. Steady-state, and time dependant. When
% doing a time dependant simulation the SmartLUSolver class should be used.
% Our simulation here is really a steady-state, but we emulate it with
% SmartLUSolver by running it for a lot of time steps.
solver = SmartLUSolver(S, es);
solver.setftoMaxmellian(); % Start with a maxwellian distribution.
solver.setnSaves(10);
output = solver.takeTimeSteps(); % Execute the simulation.

% Figure 2a
% This somewhat resembles the figures in the article.
% Note that the contour plot in the article uses units of y [p / mv_e],
% while the plotting script included here instead uses w [p / mc].
% The contour plot also includes a few additional lines which indicate
% runaway regions.
CODEPlotter.contourDistribution(output.distributions{end});

% Figure 1a
figure;
[plotF, ploty] = ParallelDistribution(output.distributions{end});
semilogy(ploty, plotF);
ylim([1e-10; 1]);
xlim([-10; 15]);

%% Figure 2
% Just plot different outputs from output.distributions.
figure;
hold on;
for dist = output.distributions
    disp(dist)
    [posF, negF] = SumLegModesAtXi1(dist{1}.f,Ny,Nxi);
    plot([-flip(S.momentumGrid.y) S.momentumGrid.y], [flip(negF); posF]);
    ylim([1e-10; 1]);
    xlim([-10; 15]);
end
hold off;
set(gca, 'YScale', 'log')

%% Figure 3
% In the article, the growthrate is calculated using 
E_space = linspace(0.04, 0.4, 10);
Z_space = [1 2 3 10];

[EE, ZZ] = meshgrid(E_space, Z_space);

EHatOverE = 0.651;
growthrates = zeros(size(EE));
for idx = 1:numel(EE)
    S.physicalParams.setE(EE(idx) * EHatOverE);
    S.physicalParams.setZ(ZZ(idx));
    solver = SolverSteadyState(S, es);
    [F, a] = solver.steadyState2();
    growthrates(idx) = a;
end

clf
hold on;
for idx = 1:length(Z_space)
    plot(E_space, growthrates(idx, :));
end
set(gca, 'YScale', 'log');
hold off;
