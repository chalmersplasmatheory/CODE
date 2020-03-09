DiagonalPlusBoundaryCondition
========================================

.. class:: DiagonalPlusBoundaryCondition

Builds the diagonal and boundary condition matrix in CODE version 1.03 CODEtimeDependent(). 
Basically does the BC and
diag matrix.
Extends ImplicitOperators but is not really a implicit operator.
Should maybe be given a new class.

Properties
-----------------
These properaties are saved todetermine if a new matrix is needed or not in next timestep.

.. attribute:: NxiUsed

.. attribute:: ddyUsed

.. attribute:: RelLimitUsed

.. attribute:: yMaxBC

Functions
---------------

.. function:: DiagonalPlusBoundaryCondition(state,eqSettings)

Creates a new instance of this class

.. function:: matrixHasChanged = generateOperatorMatrix(this,iteration)

Generates the diagonal and boundary condition from time derivative as if it was a ImplicitOperator.
