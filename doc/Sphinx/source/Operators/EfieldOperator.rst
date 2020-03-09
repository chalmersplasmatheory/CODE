EfieldOperator
===============
.. class:: EfieldOperator

Extends :class:`ImplicitOperator`. Describes the electric field term in the kinetic equation.


Functions
----------

.. function:: generateOperatorMatrix(this,iteration)
        :noindex:

generates the operator matrix for the electric field at iteration from state object, 
with Legendre mode and number of momentum points defined in state object.
Only recomputes if it relevant properties has changed.
