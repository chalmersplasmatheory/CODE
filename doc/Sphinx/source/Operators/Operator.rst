Operator
==================

.. class:: (Abstract) Operator < handle & matlab.mixin.Copyable

Abstract superclass for operators. 
Both implicit and explicit.

Properties 
-----------------------

.. attribute:: state

.. attribute:: eqSettings

.. attribute:: operatorMatrix

Sparse properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

used for sparse matrix building in CODE, not requiered for
implementation of the class, usage should be investigated in
implementations of the class

.. attribute:: predictedNNZ

.. attribute:: sparseCreatorIndex

.. attribute:: estimated_nnz

.. attribute:: sparseCreator_i

.. attribute:: sparseCreator_j

.. attribute:: sparseCreator_s

Functions
---------------------

.. function:: this = Operator(state, eqSettings, varargin)

Construct a new instance of this class.

.. function:: matrixHasChanged = generateOperatorMatrix(this,runIndex, inputArg)

generateOperatorMatrix implementation should modify operatorMatrix
to match the parameters at time index runIndex in state object 

.. function:: M = getMatrix(this)
    
The functions below are all utilities for building sparse matrices:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. function:: resetSparseCreator(this)

Resets the sparse values.

.. function:: clearSparseResidue(this)

Reduce the memory footprint by removing the vectors used for
building the sparse matrix once they are no longer needed

.. function:: addToSparse(this,i,j,s)

Adds values to the sparse matrix s at index (i,j) to matrix

.. function:: addSparseBlock(this,rowIndices, colIndices, block)

Adds a block to the sparse matrix.
rowIndices and colIndices should be vectors.
numel(rowIndices) should equal the number of rows in 'block'.
numel(colIndices) should equal the number of columns in 'block'.

.. function:: sparseMatrix = createSparse(this)

After you are done adding elements to the sparse matrix using
addToSparse() and addSparseBlock(), call this function to
finalize the matrix.
