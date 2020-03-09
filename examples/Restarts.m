%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These examples will guide you through how to restart your simulations
% It will assume that you have already created a State S, EquationSettings es and a SmartLUSolver solver.
% For these simulations, initial tMax of 10 and dt of 1 was used.
% To continue the simulation from end of takeTimeSteps and everything else is same as before, just alter tMax
% and invoke the function takeTimeSteps
S.timeGrid.settMax(200);
solver.takeTimeSteps(output);
%%
% The distributions saved will be added to the already saved cell (nothing will be overwritten).
% Moreover, all matricies are already built and therefore don't need to be rebuilt. If, for example timeGrid changes more
% drastically, timeIndex is not refering to the same time as before. Therefore use setTime function instead. setTime also can
% take 's' (seconds) or 'ms' as an additional input if the the time you want to start at is known in SI units.

S.timeGrid.setdt(0.5);
S.timeGrid.settMax(250);
solver.setTime(200,'normalized');%previous time
solver.takeTimeSteps(output);
%%
% To change the temperature if you are unsure, follow the instructions for changing the momentum grid. It is also 
% possible to just change temperature at will, in this case to from constant 10 to 15 and then continue the 
% simulation. 
T = [100, 50, 20, 10];
tT = [250, 260, 270, 280];
S.physicalParams.setinterpolationMethod('nearest');
S.physicalParams.setT(T, tT);

S.timeGrid.settMax(300)

solver.takeTimeSteps(output);

% note that if you do S.physicalParams.setT(15), it will indeed continue the simulation as in the code above. However 
% it will change the temperature in State to constant 15, which means it would seem from stat that you have simulated with T = 15
% all along. But the information of what the temperature was when the distribution was saved in the output is unchanged
%%
% To change the momentumGrid, the distributions are dependent on the already existing momentum grid, and therefore
% it is most practical to create a new state and a new solver, and then using the final output distribution to be the
% first distribution.

S2 = State(TRef,nRef);
S2.momentumGrid.setResolution('Ny',200,'yMax',50,'Nxi',8) 
S2.timeGrid.setResolution('tMax',20,'dt',1)
S2.physicalParams.setParams('E',0.2,'T',10,'n',1e19,'Z',1)
es2 = EquationSettings(S2);
% you can also mimic the settings used in another settings object by invoking the mimic function. 
% This will also copy the matricies used
%es2.mimic(es);
es2.setsourceMode(0);
es2.setcollisionOperator(0);
es2.setbremsMode(0);

solver2 = SmartLUSolver(S2,es2);
solver2.setf(output.distributions{end}, S2.momentumGrid)
solver2.takeTimeSteps();

%% it is also not necessary to use the end of the saved distribution as the new distribution, but can also be restarted from an arbitrarly 
% point in the run


S2 = State(TRef,nRef);
S2.momentumGrid.setResolution('Ny',200,'yMax',50,'Nxi',8) 
S2.timeGrid.setResolution('tMax',20,'dt',1)
S2.physicalParams.setParams('E',1,'T',10,'n',1e19,'Z',1)
es2 = EquationSettings(S2);
es2.mimic(es);
solver2 = SmartLUSolver(S2,es2);
solver2.setf(output.distributions{10}, S2.momentumGrid)
%solver2.setTime(S.distributions{10}.time); %This set so that the time steps start at time from distribution 10. 
solver2.takeTimeSteps();