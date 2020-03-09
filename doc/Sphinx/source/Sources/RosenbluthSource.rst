RosenbluthSource
================

Extends :class:`Source`.
Describes and calculates the source vectors for Rosenbluth-Putvinski avalanche source.

.. class:: RosenbluthSource

Class describing a Rosenbluth Putvinskii source. Is :attr:`sourceMode` 1 and 2 in :class:`EquationSettings`.

Properties (non-inherited)
-------------------------------

.. attribute:: tail_mask

Mask describing the tail of of the runaways. It is used in :attr:`sourceMode` 2 from :class:`EquationSettings`.

.. attribute:: nr_mask

Mask defining the runaways, i.e. when an electron is considered a runaway. This is done so all with energy bigger twice the critical energy is considered to be a runaway (if not wMinForSource is set by user).

Functions
------------

.. function:: getSourceVec(this, f, iteration)

Returns the source vector for this iteration and also updates the :attr:`sourceVector`.

.. function:: DetermineFastParticleRegion

Internal function whihc determines the tail mask.
