%Checks if matricies are the same and if timestep vector along with
%momentum grid is the same

T = 10^(rand()*3);
n = 10^(rand()*1.4)*1e19;
Z = rand*4+1;
E = rand*2;
B = 0;%(rand()>0.4)*(rand*6+2) + (rand<0.02)*200;
% species                      = struct();
% species.Z0NetCharge          = [1;1]; % Hydrogen and Argon
% species.ZAtomicNumber        = [1;18];
% species.nj      = zeros(length(species.Z0NetCharge), 1);
% species.nj(1,:) = 1e19;
% species.nj(2,:) = 1e18;
% species.times   = 0;
Nxi = 50;
Ny = 2000;
yMax = 150;
dt = rand*2+0.1;
tMax = rand*50+50;
TRef = T;
nRef = n;
timeStepMode = 0;
%timeStepMode = floor(rand()*4);
if timeStepMode == 1
    timeStepMode = 0;
end
logGridScaling = rand+1;
logGridSubSteps = round(rand*5+2);
logGridMaxStep = rand*50+200;
collisionOperator = 0;%floor(rand*5);
enforceDensityConservation = 1;%round(rand);
bremsMode = 0;%floor(rand*5);
NyInterp = round(rand*50+200);
useScreening = 0;%round(rand);
useInelastic = 0;
sourceMode = 0;%floor(rand*6);
fastParticleDefinition = 0;%floor(rand*5);
tailThreshold = rand*4+3;
yCutSource = rand*4+3;
yGridMode = 6;%floor(rand*7);
if yGridMode == 5 || yGridMode == 2
    yGridMode = 6;
end
gridParameter = (yGridMode == 6)*(rand*2+1) + (yGridMode == 0 || yGridMode == 4)*rand*0.1;
gridStepWidth = (yGridMode == 6)*(rand*2+4) + (yGridMode ~= 6 )*((0.5-rand)*2/50+1/50);
gridStepPosition = (yGridMode == 6)*(rand*0.9+0.05) + (yGridMode ~= 6 )*((0.5-rand)*1/2+2);
yMaxBoundaryCondition = 0;%floor(rand*4); %floor(rand*5);
artificialDissipationStrength = 1e-2;%artificialDissipationStrength = rand*0.1;
artificialDissipationWidth = 1e-1;%artificialDissipationWidth = rand*20+1;
nStepsToReturn = tMax/dt +10;%floor(rand*30+20);
useFullSynchOp = 1; %useFullSynchOp = round(rand);
useNonRelativisticCollOp = 0;%round(rand);
ipm = {'nearest','linear','spline'};
interpolationMethod = 'linear';ipm{floor(rand*2+1)};
runawayRegionMode = 0;%round(rand);
nPointsXiInt = 30+round(rand*40);
autoInitial = 0;

% Time dependence
% tT = 0:dt:tMax;
% T = T * exp(-tT/tMax);
% n = n * exp(-tT/tMax);
% E = E * exp(-tT/tMax);
% Z = Z * exp(-tT/tMax);


%%%%%%%%%%%%%%%%%%%%%%%% Old CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o = CODESettings();
o.keepTrackOfExternalRunaways = 0;
o.T = T;
o.n = n;
o.Z = Z;
o.E = E;
o.B = B;
%%%%%%%%%%%%% Time dep %%%%%%%%%%%%%5
% o.tT = tT;
% o.tn = tT;
% o.tZ = tT;
% o.tE = tT;

%%%%%%%%%%%% Species %%%%%%%%%%%%%%%%
%o.species = species;

o.Nxi = Nxi;
o.Ny = Ny;
o.yMax = yMax;
o.dt = dt;
o.tMax = tMax;
o.TRef = TRef;
o.nRef = nRef;
o.timeUnit = 'normalized';
o.timeStepMode = timeStepMode;
o.logGridScaling = logGridScaling;
o.logGridSubSteps = logGridSubSteps;
o.logGridMaxStep = logGridMaxStep;
o.collisionOperator = collisionOperator;
o.enforceDensityConservation = enforceDensityConservation;
o.bremsMode = bremsMode;
o.NyInterp = NyInterp;
o.useScreening = useScreening;
o.useInelastic = useInelastic;
o.sourceMode = sourceMode;
o.fastParticleDefinition = fastParticleDefinition;
o.tailThreshold = tailThreshold;
o.yCutSource = yCutSource;
o.yGridMode = yGridMode;
o.gridParameter = gridParameter;
o.gridStepWidth = gridStepWidth;
o.gridStepPosition = gridStepPosition;
o.yMaxBoundaryCondition = yMaxBoundaryCondition;
o.artificialDissipationStrength = artificialDissipationStrength;
o.artificialDissipationWidth = artificialDissipationWidth;
o.nStepsToReturn = nStepsToReturn;
o.useFullSynchOp = useFullSynchOp;
o.useNonRelativisticCollOp = useNonRelativisticCollOp;
o.interpolationMethod = interpolationMethod;
o.runawayRegionMode = runawayRegionMode;
o.nPointsXiInt = nPointsXiInt;
o.silentMode = true;
o.reduceMemoryConsumption = false;
o.autoInitialGrid = autoInitial;

%res = CODE_timeDependent(o);
[CODE_Operator, CODE_source_vec, CODE_heat_part_preserv_vec, CODE_timesteps, CODE_momentumGrid] = CODE_timeDependent_returnMatricies(o); %this is basically a code that returns its operators and rhs after its first iteration


%%%%%%%%%%%%%%%%%%%%%%%% CODE Structured %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = State(TRef,nRef);
%%%%%%%%%%%% time dep or not %%%%%%%%%%%%%%%%%
S.physicalParams.setT(T);%,tT);
S.physicalParams.setn(n);%,tT);
S.physicalParams.setZ(Z);%,tT);
S.physicalParams.setE(E);%,tT);
S.physicalParams.setB(B);
S.physicalParams.setinterpolationMethod(interpolationMethod);
%%%%%%%%%%% species %%%%%%%%%%%%%%%%
%S.physicalParams.setspecies(species);
S.timeGrid.setdt(dt);
S.timeGrid.settMax(tMax);
S.timeGrid.settimeStepMode(timeStepMode);
S.timeGrid.setlogGridScaling(logGridScaling);
S.timeGrid.setlogGridSubSteps(logGridSubSteps);
S.timeGrid.setlogGridMaxStep(logGridMaxStep);
S.momentumGrid.setNxi(Nxi);
S.momentumGrid.setNy(Ny);
S.momentumGrid.setyMax(yMax);
S.momentumGrid.setyGridMode(yGridMode);
S.momentumGrid.setgridParameter(gridParameter);
S.momentumGrid.setgridStepWidth(gridStepWidth);
S.momentumGrid.setgridStepPosition(gridStepPosition);
S.momentumGrid.setyMaxBoundaryCondition(yMaxBoundaryCondition);
S.momentumGrid.setartificialDissipationStrength(artificialDissipationStrength);
S.momentumGrid.setartificialDissipationWidth(artificialDissipationWidth);
if autoInitial
    S.autoInitGrid(useScreening,useInelastic)
end

es = EquationSettings(S);
es.setcollisionOperator(collisionOperator);
es.setenforceDensityConservation(enforceDensityConservation);
es.setbremsMode(bremsMode);
es.setuseScreening(useScreening);
es.setsourceMode(sourceMode);
es.setfastParticleDefinition(fastParticleDefinition);
es.settailThreshold(tailThreshold);
es.setyCutSource(yCutSource);
es.setNyInterp(NyInterp);
es.setuseFullSynchOp(useFullSynchOp);
es.setuseNonRelativisticCollOp(useNonRelativisticCollOp);
es.setnPointsXiInt(nPointsXiInt);
es.setrunawayRegionMode(runawayRegionMode);
es.setuseInelastic(useInelastic)

solver = SmartLUSolver(S,es);
solver.setnSaves(nStepsToReturn);
solver.setftoMaxmellian()

N_CODE_Operator = sparse(Nxi*Ny,Nxi*Ny);
for IO = solver.implicitOperators
    IO{1}.generateOperatorMatrix(1);
    if isa(IO{1},'DiagonalPlusBoundaryCondition')
        %N_CODE_Operator = N_CODE_Operator + IO{1}.getMatrix();
    else
        N_CODE_Operator = N_CODE_Operator + IO{1}.getMatrix();
    end
end

N_CODE_source_vec = zeros(Nxi*Ny,1);
for EO = solver.explicitOperators
    EO{1}.generateOperatorMatrix(1);
    N_CODE_source_vec = N_CODE_source_vec + solver.state.timeGrid.dts(1)*EO{1}.getMatrix*S.fattimeIndex;
end
for source = solver.sources
    N_CODE_source_vec = N_CODE_source_vec + solver.state.timeGrid.dts(1)*source{1}.getSourceVec(S.fattimeIndex,1);
end
N_CODE_heat_part_preserv_vec = solver.enforceParticleAndHeat();
N_CODE_timesteps = S.timeGrid.timesteps;
N_CODE_momentumGrid = S.momentumGrid.y;

disp('Max relative diff Operator')
diff = N_CODE_Operator+CODE_Operator;
idx = find(CODE_Operator);
disp(max(max(abs(diff(idx)./CODE_Operator(idx)))))
disp('Max diff source_vec')
disp(max(max(abs(N_CODE_source_vec-CODE_source_vec))))
disp('Max diff heat_part_preserv_vec')
disp(max(max(abs(N_CODE_heat_part_preserv_vec-CODE_heat_part_preserv_vec))))
disp('Max diff timesteps')
disp(max(max(abs(N_CODE_timesteps-CODE_timesteps))))
disp('Max diff momentumGrid')
disp(max(max(abs(N_CODE_momentumGrid-CODE_momentumGrid))))
