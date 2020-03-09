classdef Source < handle
    %SOURCE Abstract class defining a Source for a CODE solver. Contains
    %grid, plasma and source properties. Is used when operators are hard
    %or infeasible to write as a matrix operator.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%        Properties                %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%
    %  General Properties:
    % %%%%%%%%%%%%%%%%%%%%%%%%
    %
    % state - State object containing grid and physical quantaties
    % sourceVector - vector with the source

    properties
        state
        eqSettings
        sourceVector
    end

    properties
       VERSION = 1.0;
    end

    methods
        function this = Source(state, eqSettings)
            %SOURCE Constructor defining a source for a grid and plasma
            %object
            this.state = state;
            this.eqSettings = eqSettings;
        end
    end


    methods (Abstract)
        % Returns the source vector for the source. 
        getSourceVec(this,f,iteration)
    end
    
end

