Sources
==============

.. class:: Source

Abstract class defining a avalanche Source. Can for example be the avalanche source Rosenbluth and Putvinskii derived.

Properties
------------------------

.. attribute:: state

:class:`State` object for knowing the physical params use for knowing the physical params usedd

.. attribute:: eqSettings

:class:`EquationSettings` object used for knowing what settings is used

.. attribute:: sourceVector

sourceVector - vector with the source, calculated each step with getSourceVec.

Functions
---------------

.. function:: getSourceVec(this,f,iteration)

Returns the source vector at iteration from f.
Should recalc and update :attr:`sourceVector`.
