%% Figure 1a
n = 1e20; % [m^-3]
T = 10; % [eV]
B = 5;
Z = 1;

%species.times = 1; % What 
species.times = 0;
species.nj = [n; n];
species.ZAtomicNumber = [1; 18];
species.Z0NetCharge = [1; 1];

% Calculate E_D:
eps_0 = 8.854187e-12;                   % As/Vm
e = 1.60218e-19;                        % As
eV = 1.60218e-19;                       % J

ne = sum(species.nj .* species.Z0NetCharge);
log_lambda = 14.9 - 0.5 * log(ne / 1e20) + log(T / 1e3);
E_D = ne * e^3 * log_lambda / (4 * pi * eps_0^2 * T) / eV;

Espace = linspace(0, 0.1, 10) * E_D;

% Since we will not change T och n during the simulation, we set the
% reference parameters to the same thing. The reference parameters
% determines how the grid is normalized.
S = State(T, ne);
S.physicalParams.setPhysicalParameters(T, ne, 1, Espace(1));
S.physicalParams.setspecies(species);

% The same resolution paramters are used as in the article.
yMax = 20;
Ny = 300;
Nxi = 20;
S.momentumGrid.setResolution('yMax', yMax, 'Ny', Ny, 'Nxi', Nxi);

% CODE is very efficient when the physical paramters aren't changed during
% the simulation. Which is why we can set tMax very high here.
tMax = 500;
S.timeGrid.setResolution('tMax',tMax,'dt',1);

% There are two solvers in CODE. Steady-state, and time dependant. When
% doing a time dependant simulation the SmartLUSolver class should be used.
% Our simulation here is really a steady-state, but we emulate it with
% SmartLUSolver by running it for a lot of time steps.
screening = [0 0 1];
inelastic = [0 0 1];
edlnlambda = [0 1 1];

growthrates = zeros(3, numel(Espace));
for l=1:3
%growthrates = zeros(size(Espace));
for k = 1:numel(Espace)
    E=Espace(k);
es = EquationSettings(S);
es.setuseScreening(screening(l));
es.setuseInelastic(inelastic(l));
es.setuseEnergyDependentLnLambdaScreening(edlnlambda(l));

S.physicalParams.setE(E);
solver = SmartLUSolver(S, es);
solver.setftoMaxmellian(); % Start with a maxwellian distribution.
solver.setnSaves(1);
output = solver.takeTimeSteps(); % Execute the simulation.

growthrates(l, k) = GetFlowThroughGridBoundary(output.distributions{end}, yMax - 4) * S.physicalParams.nuees(end);
end
end

plot(Espace / E_D, growthrates(1, :));
hold on;
plot(Espace / E_D, growthrates(2, :));
plot(Espace / E_D, growthrates(3, :));
hold off;
set(gca, 'YScale', 'log');
ylim([1e-7 1e6]);
