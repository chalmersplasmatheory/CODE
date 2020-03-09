ImplicitOperator
=================
.. class:: ImplicitOperator

This is and abstract class handling implicit operators. These can ofcourse be used as explicit operators by adding a minus sign to the matrix after generation of the matrix.

Properties
------------

.. attribute:: state

:class:`State` object for creating the matricies 

.. attribute:: eqSettings

:class:`EquationSettings` object, containing what settings to use for the run.

.. attribute:: operatorMatrix

the operatorMatrix. Is calculated with the function call :func:`generateOperatorMatrix`.

Sparse properties:
^^^^^^^^^^^^^^^^^^^^^^^

used for sparse matrix building in several implimations of this class. It is not requiered for
implementation though.

.. attribute:: predictedNNZ 0;
.. attribute:: sparseCreatorIndex=1;
.. attribute:: estimated_nnz = 0;
.. attribute:: isparseCreator_i=0;
.. attribute:: sparseCreator_j=0;
.. attribute:: sparseCreator_s=0;


Functions
------------

.. function:: generateOperatorMatrix(this,runIndex, inputArg)

generateOperatorMatrix implementation should modify operatorMatrix to be the operator matrix at runIndex.
Implementation should also return 1 (true) if matrix changed from previously saved matrix and 0 otherwise.
