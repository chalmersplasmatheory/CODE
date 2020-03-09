ImprovedChiuHarveySource
==========================
.. class:: ImprovedChiuHarveySource

Calculates pitch angle dependent Chiu Harvey like source and can be switched to have a particle conserving sink term to consverve number of particles and momentum (arissen from the source term).
Extends :class:`ExplicitOperator`.
Is :attr:`sourceMode` 4,5 in :class:`EquationSettings`.

Functions
----------------

.. function:: generateOperatorMatrix(this, iteration)

Calculates the explicit source term at the iteration as the time index used.
Uses paramters from the :class:`State` object.
Returns true if matrix has changed and 0 otherwise.

.. function:: GenerateKnockOnMatrix(this,iteration)

Generates the KnockOn terms for the ChiuHarveySource with pitch angle dependencies.

.. function:: BuildInterpolationMatrix(~,x,xp)

Help function used only in the class.
