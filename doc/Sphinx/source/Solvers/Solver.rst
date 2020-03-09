Solver
===============

.. class:: Solver

Abstract solver class structuring a solver class to be used in CODE

Properties
-------------

.. attribute:: state

State object containing physical information about the
plasma such as temperature, Efield, Bfield and about what momentum grid
is used. This is a handle class shared by all operator
objects

.. attribute:: implicitOp

cell containing all implicit operators used in the run
all objects in cell should extend ImplicitOperator class

.. attribute:: explicitOp

cell containing all explicit operators used in the run
all objects in cell should extend ExplicitOperator class

.. attribute:: sources

cell containing all sources for the run, care these are
not decided how they should be structured inside
the code yet

.. attribute:: eqSettings

EquationSettings containing what settings for the
equation to use

Functions
----------------

.. function:: this = Solver(state, eqSettings)

Create a new isntance of this class.


.. function:: updateOperators(this)

Updates the operators in this class to match those of the :attr:`eqSettings` (:class:`EquationSettings` object).

.. function:: takeTimeSteps(this)

Abstract function, should step in time somehow.
