ExplicitOperator
======================

.. class:: (Abstract) ExplicitOperator < handle & Operator

EPLICITOPERATOR an abstract superclass for explicit operators used in CODE.
Defining abstract functions for updating its operation vector and non abstract utility functions for creating sparse matricies.

Properties
----------------

VERSION = 1.0

Functions
----------------

function this = ExplicitOperator(state, eqSettings, varargin)

Construct a new instance of this class

