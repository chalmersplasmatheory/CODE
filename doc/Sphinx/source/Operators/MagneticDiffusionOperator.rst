MagneticDiffusionOperator
================================

.. class:: MagneticDiffusionOperator < ImplicitOperator

Calculates local transport sink.
Extends :class:`ImplicitOperator`.

Properties
--------------------

.. attribute:: deltaBoverB

Size of radial magnetic field fluctuations

.. attribute:: safetyFactor

Safety factor (q) needed for diffusion coefficient estimate

.. attribute:: majorRadius

Major radius (in meters) for Rechester-Rosenbluth diffusion

.. attribute:: radialScaleLength

Length scale for radial variation of runaway density (in meters)

.. attribute:: ionMassNumber

Ion mass number (in proton masses)

Functions
====================

function this = MagneticDiffusionOperator(state, eqSettings, varargin)

Creates a new instance of this class.

function matrixHasChanged = generateOperatorMatrix(this,iteration)
