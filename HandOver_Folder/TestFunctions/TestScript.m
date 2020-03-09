%% 
months = {{'Jan'}, {'Feb'}, {'Mar'}, {'Apr'}, {'May'}, {'Jun'}, {'Jul'}, {'Aug'}, {'Sep'}, {'Oct'}, {'Nov'}, {'Dec'}};
% successfull
results = struct('filename', '','currDens',[], 'Energy',[], 'Distribution',[], 'MomentumGrid',[],'Time',[]);
for year = 2019
    for month = {months{10}}
        month = month{1,1};
        month = char(month);
        for day = 9
            for i = 1
                filename = ['TestResults/Test_' num2str(day,'%02.f') '-' month '-' num2str(year) '_' num2str(i) '.mat'];
                if isfile(filename)
                    load(filename)
                    resCODE.times(end)-resCODENew.distributions{end}.time
                    %resCODE.settings.dt
%                     if (length(resCODENew.averageEnergyMeV)-resCODE.settings.nStepsToReturn) ~= 0
%                         return
%                     end
                    norm((resCODE.f(:,end)-resCODENew.distributions{end}.f))/length(resCODE.f(:,end))
                    results(end+1) = struct('filename',filename,...
                        'currDens', 1e-60 > abs((resCODE.currentDensity(end)-resCODENew.currentDensity{end})),...
                        'Energy',  1e-6 > abs((resCODE.averageEnergyMeV(end)-resCODENew.averageEnergyMeV{end})), ...
                        'Distribution', 1e-60 > norm((resCODE.f(:,end)-resCODENew.distributions{end}.f)), ...
                        'MomentumGrid', 1e-40 > norm((resCODE.y - resCODENew.distributions{end}.momentumGrid.y)),...
                        'Time', resCODE.times(end)-resCODENew.distributions{end}.time);
                end
            end
        end
    end
end
pause(1)
disp(length(find([results.currDens]))/length([results.currDens]))
disp(length(find([results.Energy]))/length([results.currDens]))
disp(length(find([results.Distribution]))/length([results.currDens]))
disp(length(find([results.MomentumGrid]))/length([results.currDens]))
disp([results.Time])
CODEDIST = [resCODENew.distributions{:}];
resCODE.times - [CODEDIST.time]
%% Error

for year = 2019
    for month = months{9}
        month = char(month);
        for day = 3
            for i = 1:30
                filename = ['TestResults/ErrorCODE_' num2str(day,'%02.f') '-' month '-' num2str(year) '_' num2str(i) '.mat'];
                if isfile(filename)
                    load(filename)
                    if exist('e2','var')
                        e2.getReport
                    end
                    if exist('e','var')
                        e.getReport
                    end
                    clear 'e2' 'e'
                end
            end
        end
    end
end

%%
notSame = {results(find(~[results.Distribution])).filename};

load(char(notSame{2}))
%%
norm(resCODE.f(:,i) - resCODENew.distributions{i}.f)
resCODE.f(1:10,i)
resCODENew.distributions{i}.f(1:10)
    

%%

if resCODE.settings.T ~= resCODENew.T{1}
disp('T');
end
if resCODE.settings.n ~= resCODENew.n{1}
disp('n');
end
if resCODE.settings.E ~= resCODENew.E{1}
disp('E');
end
 if resCODE.settings.Z ~= resCODENew.Z{1}
disp('Z');
end
 if resCODE.settings.B ~= resCODENew.B{1}
disp('B');
end
% if resCODE.neTotalOverneFree ~= resCODENew.neTotalOverneFree{1}
% disp('neTotalOverneFree');
% end
 if resCODE.nuee ~= resCODENew.nuees{1}
disp('nuees');
end
 if resCODE.delta ~= resCODENew.deltas{1}
disp('deltas');
end
 if resCODE.EOverEc ~= resCODENew.EOverEc{1}
disp('EOverEc');
end
if resCODE.EOverED ~= resCODENew.EOverED{1}
disp('EOverEd');
end
 if resCODE.EHat ~= resCODENew.EHats{1}
disp('EHats');
end
 if resCODE.BHat ~= resCODENew.BHatRef{1}/sqrt(resCODENew.nueeBars{1})
disp('BHatRef');
end
 if resCODE.nuee/resCODE.nueeRef ~= resCODENew.nueeBars{1}
disp('nueeBars');
end
if resCODE.nBars ~= resCODENew.nBars{1}
disp('nBars');
end
if resCODE.veBars ~= resCODENew.veBars{1}
disp('veBars');
end
% if resCODE.lnLambdas ~= resCODENew.lnLambdas{1}
% disp('lnLambdas');
% end

