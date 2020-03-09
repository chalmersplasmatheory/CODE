SynchrotronLoss
=====================

.. class:: SynchrotronLoss

Class used to calculate the syncrotron radiation from disitribution.
Extends the implicit operator.

Properties
------------------

Properties used for the given operator
*************************************************

.. attribute:: nueeBarUsed

What collision frequency is used for creating matrix.

.. attribute:: BHatUsed

What normalized magnetic field is used for creating matrix.

.. attribute:: NxiUsed

What number of legandre mode is used for creating matrix.

.. attribute:: yUsed

What momentum grid is used for creating matrix.

.. attribute:: deltaRefUsed

Which normalization factor is used for creating matrix.

.. attribute:: gammaUsed

What gamma is used for creating matrix.

.. attribute:: yMaxBCUsed

What boundary condition is used for creating matrix. (Relevant for where the matrix should start)

.. attribute:: useFullSynchOpUsedA

If full synchrotron operator is used for creating matrix.

Functions
-----------------

.. function:: this = SynchrotronLoss(state,eqSettings)

Creates a new instance of this class.

.. function:: matrixHasChanged = generateOperatorMatrix(this,iteration)
