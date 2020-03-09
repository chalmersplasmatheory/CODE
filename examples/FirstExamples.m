%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Following will guide you through how to write a simple runscript
%
% The code is seperated into three main components: Solvers, Operators and State (and a fourth, namely the output)
% The Solver is a class solving your specific setup. The 'state' of the plasma, i.e.
% its temperature, the electric field, magnetic field and more; together with what timepoints
% is of intrest; at which resolution your momentum grid will be at; the normalization for numerical
% stability used. 
%
% Information about how the equation is built up is contained inside the Operators, and
% is collect in an instance of the EquationSettings class. This class then contains all 
% relevant operators for your scenario.
%
% To start, we first need a State, in which the plasma is in. This in turn needs to know how it is normalized.
TRef = 10;
nRef = 1e19;
S = State(TRef,nRef);
% Then we can set the relevant properties of the State by calling different set functions.
% The State S now contains 4 1x1 cells of objects and one cell where distribution functions will be stored

%% Momentum
% To set things related to the momentum, we use the S.momentumGrid. This contains set functions with names
% 'set[Param]' where [Param] is replaced with what parameter you want to use. It also contains the function
% setResolution which takes input as ('name1',value1,'name2',value2 ...) and sets the respective values.
% The momentum grid is now initiated
S.momentumGrid.setResolution('yMax',200,'Ny',600,'Nxi',40,'yGridMode',4)

%% TimeGrid
% The time grid works exactly as the momentum grid, with the exception that it is accessed from S.timeGrid
% The unit supplied will always be in normalied time, that is normalized to the reference collision time.
% This is located in S.reference, and it is therefore as easy to convert your time from seconds to normalized 
% time by the simple call
tMax = 1e-3*S.reference.nueeRef;%simulation will proceed until 1 ms
S.timeGrid.setResolution('tMax',tMax,'dt',1)

%% It is also possible to create a own time vector to use and then invoke settimesteps. This takes takes the 
% timesteps vector and what unit it is supplied in. Possible are 's' (seconds) 'ms' (milliseconds) and 
% 'normalized' (normalized).
% This can be useful for example when temperature or density change so collision times changes and therefore 
% the dynamics and how much resolution is required in time during different times. By not restarting the matricies
% will not be rebuilt, only the inversion will be redone when dt changes, therefore speeding up compared to restarting
% This is sofar not compatibly with all restarts possibilities (but is the
% only thing not compatible sofar)
timesteps = [0:0.1:10, 10:1:20, 20:5:100];
S.timeGrid.settimesteps(timesteps, 'normalized')

%% If you now want to actually used an built in timestep mode we need to
% change timstepmode
% lets use a homogeneous timegrid
S.timeGrid.settimeStepMode(0)
S.timeGrid.setdt(0.1)
S.timeGrid.settMax(100)

%% Physical Quantities
% Physical Quantities is located in S.physicalParams (params for parameters). To set these, either call them one by one
% with 'set[Param]' where [Param] is replaced with name of param (ex setT). Another option is to call the function setParams
% which takes input as ('name1',value1,'name2',value2 ...). A third option for old CODE users is the setPhysicalParameters function
% which takes input in one of the following formats:
% Usage:
% S.physicalParams.setPhysicalParameters(T,n,Z,E)
% S.physicalParams.setPhysicalParameters(T,n,Z,E,B)
% S.physicalParams.setPhysicalParameters(T,n,Z,E,B,tT,tn,tZ,tE)
% The t* vectors are, just as in the case of TimeGrid in normalized time units.
% Other units are: T in eV, n in m^-3, Z in negativly elementary charges, E in V/m and B in Teslas.

S.physicalParams.setPhysicalParameters(20,1e19,1,0.5)

% if the temperature is time dependent, then it can be set from setT
S.physicalParams.setT(1000:-18:100,0:0.1:5)

%% If we now realize that something is wrong, for example the normalization (or just want to change it)
% it is now possible, and everything in state will update accordingly
% All time vectors will be same in absolute time, but will ofcourse change in normalized units as these have changed.
S.reference.updateReferenceVals(20,1e19)


%% Operators
% It is now time to create the operators used for the run, and what settings are used. To know what each setting do, read the documentation.
% If something is missing in the documentation, check the old CODE documentation
% The operators will be created according to what settings used. To be able to use the settings we need a EquationSettings object.
es = EquationSettings(S);

% To altare settings, we use the set functions named 'set[SettingName]'. 
es.setsourceMode(0);
es.setcollisionOperator(0);
es.setbremsMode(0);

%% Solver
% We are now ready to solve the equation, and is simply done by creating a solver and invoking the takeTimeSteps function
% To introduce what initial condition on the distribution we need to set the distribution at current time index in solver.
% This can either be done with a Disitribution object and an optional MomentumGrid object to interpolate the distribution to
% by invoking the function setf(Distribution, MomentumGrid(optional)),or by setting it to a maxwellian distribution by simply 
% invoking setftoMaxwellian();
% How many distributions are saved is altered by the nStepsToSave variable and is set with setnStepsToSave;
solver = SmartLUSolver(S,es);
solver.setftoMaxmellian();
solver.setnSaves(50);
output = solver.takeTimeSteps();

% The program has now been run and 50 steps are saved. These are saved in your state variable in a cell array called distributions and 
% are saved as Distribution objects.

% Plotting
% to plot the result, it is as simple as using the static class CODEPlotter. Just pass one of the distribution objects to one of its plotting function.
% To plot many distribution functions I would recomend using cellfun and making use of the figureOffset in the plotting functions
CODEPlotter.contourDistribution(output.distributions{end})
CODEPlotter.Plotp2f(figure(),output.distributions{end});%Will plot the last disitrubition 