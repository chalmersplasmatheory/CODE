classdef (Abstract) ExplicitOperator < handle & Operator
    %EPLICITOPERATOR an abstract superclass for explicit operators used in
    %CODE. Defining abstract functions for updating its operation vector
    %and non abstract utility functions for creating sparse matricies.

    properties (Constant)
       VERSION = 1.0
    end
    
    methods
        function this = ExplicitOperator(state, eqSettings, varargin)
            this@Operator(state,eqSettings)
        end
    end
end

