classdef (Abstract) ImplicitOperator < handle & Operator
    %IMPLCICITOPERATOR an abstract superclass for implicit operators used in
    %CODE. Defining abstract functions for updating its operator vector
    %and non abstract utility functions for creating sparse matricies.

    properties (Constant)
       VERSION = 1.0;
    end
    
    methods
        
        function this = ImplicitOperator(state, eqSettings, varargin)
            this@Operator(state, eqSettings, varargin)
        end
        
    end
    
end