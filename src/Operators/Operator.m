classdef (Abstract) Operator < handle & matlab.mixin.Copyable
    %OPERATOR operator class
    
    properties 
        state
        eqSettings
        operatorMatrix
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Sparse properties:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   used for sparse matrix building in CODE, not requiered for
    %   implementation of the class, usage should be investigated in
    %   implementations of the class
    %
    properties
        predictedNNZ = 0;
        sparseCreatorIndex=1;
        estimated_nnz = 0;
        sparseCreator_i=0;
        sparseCreator_j=0;
        sparseCreator_s=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%        Constructor         %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        % create a new Operator with the plasma and grid used in
        % Solver: NOTE: grid and plasma are handle classes, dont change
        % params there! Whenever an implementation is made, make sure you
        % call Operator with [object]@Operator([your state], [your equationSettings])
        function this = Operator(state, eqSettings, varargin)
            if isa(state,'State') && state.VERSION == this.VERSION
                this.state = state;
            else
                error(['The constructor only accepts Grid object of version ' num2str(this.VERSION) '.']);
            end
            if isa(eqSettings, 'EquationSettings') && eqSettings.VERSION == this.VERSION
                this.eqSettings = eqSettings;
            else
                error(['The constructor only accepts EquationSettings object of version ' num2str(this.VERSION) '.']);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%     matrix building        %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Abstract)
        %GENERATEOPERATORMATRIX - generate one matrix by size of grid as
        %the operator
        matrixHasChanged = generateOperatorMatrix(this,runIndex, inputArg)
            % generateOperatorMatrix implementation should modify operatorMatrix
            % to Grid object and PhysicalParams object, at temprature plasma.T(runIndex)
            % where plasma.T(runIndex) will be replaced with
            % plasma.*(runIndex) for every other parameter in plasma
    end
    
    methods
        function M = getMatrix(this)
            M = this.operatorMatrix;
        end
    end
    %%%%% sparse methods used for all CollisionalOperators in current
    %%%%% CodeVersion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The functions below are all utilities for building sparse matrices:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function resetSparseCreator(this)
            this.sparseCreatorIndex=1;
            this.estimated_nnz = floor(this.predictedNNZ);
            this.sparseCreator_i=zeros(this.estimated_nnz,1);
            this.sparseCreator_j=zeros(this.estimated_nnz,1);
            this.sparseCreator_s=zeros(this.estimated_nnz,1);
        end
        function clearSparseResidue(this)
            %Reduce the memory footprint by removing the vectors used for
            %building the sparse matrix once they are no longer needed

            %A: cannot clear properties in matlab, thus changing them to
            %negative values will create error index out of bounds which
            %means the code wont run with unknowing errors
            this.sparseCreator_i = -1;
            this.sparseCreator_j = -1;
            this.sparseCreator_s = -1;
        end
        function addToSparse(this,i,j,s)
            % Adds values to the sparse matrix.
            n=numel(i);
            if n ~= numel(j)
                error('Error A');
            end
            if n ~= numel(s)
                error('Error B');
            end
            if any(i<1)
                error('Error Q: i<1');
            end
            if any(j<1)
                error('Error Q: j<1');
            end
            this.sparseCreator_i(this.sparseCreatorIndex:(this.sparseCreatorIndex+n-1)) = i;
            this.sparseCreator_j(this.sparseCreatorIndex:(this.sparseCreatorIndex+n-1)) = j;
            this.sparseCreator_s(this.sparseCreatorIndex:(this.sparseCreatorIndex+n-1)) = s;
            this.sparseCreatorIndex = this.sparseCreatorIndex+n;
            if this.sparseCreatorIndex > this.estimated_nnz
                %fprintf('Warning! estimated_nnz is too small. Increase predictedFillFactor.\n')
                warning('Estimated_nnz is too small. Increase predictedFillFactor.')%this is for debugging GO+CODE, but usually it means something has gone wrong
            end
        end
        function addSparseBlock(this,rowIndices, colIndices, block)
            % Adds a block to the sparse matrix.
            % rowIndices and colIndices should be vectors.
            % numel(rowIndices) should equal the number of rows in 'block'.
            % numel(colIndices) should equal the number of columns in 'block'.
            s=size(block);
            if (s(1) ~= numel(rowIndices)) || (s(2) ~= numel(colIndices))
                error('Error in addSparseBlock!  Input sizes are not consistent.')
            end
            [rows, cols, values] = find(block);
            this.addToSparse(rowIndices(rows),colIndices(cols),values)
        end
        function sparseMatrix = createSparse(this)
            % After you are done adding elements to the sparse matrix using
            % addToSparse() and addSparseBlock(), call this function to
            % finalize the matrix.
            sparseMatrix = sparse(this.sparseCreator_i(1:(this.sparseCreatorIndex-1)), this.sparseCreator_j(1:(this.sparseCreatorIndex-1)), this.sparseCreator_s(1:(this.sparseCreatorIndex-1)), this.state.momentumGrid.matrixSize, this.state.momentumGrid.matrixSize);
            this.resetSparseCreator()
            %         warning('The sparse matrix creator is not automatically reset');
        end
    end
end

