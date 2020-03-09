.. CODE2 documentation master file, created by
   sphinx-quickstart on Wed Feb  6 11:51:19 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CODE documentation
=================================
Code documentation, its uses and internal structure

Classes
-------------
Below follows the different classes in CODE.

.. toctree::
   :hidden:

   self


Sources
%%%%%%%%%%%%%%%%%%

.. toctree::
    :maxdepth: 1 

    Sources/Source
    Sources/RosenbluthSource
    Sources/ChiuHarveySource

Operators
%%%%%%%%%%%

.. toctree::
    :maxdepth: 1

    Operators/Operator
    Operators/EquationSettings
    Operators/DiagonalPlusBoundaryCondition

Implicit Operators
%%%%%%%%%%%%%%%%%%%%%%%

.. toctree::
    :maxdepth: 1

    Operators/ImplicitOperator
    Operators/EfieldOperator
    Operators/CollisionOperator
    Operators/SynchrotronLoss
    Operators/MagneticDiffusionOperator
   
Explicit Operators
%%%%%%%%%%%%%%%%%%%%%%%

.. toctree::
    :maxdepth: 1

    Operators/ExplicitOperator
    Operators/ImprovedChiuHarveySource
    Operators/Bremsstrahlung

Quantities
%%%%%%%%%%%%%%%%%%%

.. toctree::
    :maxdepth: 1

    State/State
    State/Reference
    State/PhysicalParams
    State/MomentumGrid
    State/TimeGrid

Simulation
%%%%%%%%%%%%%%%%

.. toctree::
    :maxdepth: 1

    Solvers/Solver
    Solvers/SmartLUSolver

Simulation output
%%%%%%%%%%%%%%%%%%%%%%%%%%

.. toctree::
    :maxdepth: 1

    Output/Output
    Output/Distribution

Help functions
--------------------
Below are different help functions located in CODElib folder documented.

.. toctree::
    :maxdepth: 1

    CodeStructure/HelpFuncs

Code internal structure
-------------------------
.. toctree::
    CodeStructure/variableStructure

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
