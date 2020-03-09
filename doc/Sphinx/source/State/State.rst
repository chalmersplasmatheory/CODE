State
===========

.. class:: State 

Class gathering a collection of :class:`PhysicalParams`, :class:`MomentumGrid`, :class:`TimeGrid` and :class:`Reference` objects which point to each other according to their documentation.
Can be seen as the total state and resolution at which the equation is solved for.
Use this class to initiate the :class:`PhysicalParams`, :class:`MomentumGrid`, :class:`TimeGrid` and :class:`Reference`.

properties 
---------------

.. attribute:: VERSION 

Version of class

Dependencies
%%%%%%%%%%%%%%%%%

.. attribute:: physicalParams

:class:`PhysicalParams` object shared with timeGrid and reference

.. attribute:: reference

:class:`Reference` object shared with all other objects in this class

.. attribute:: momentumGrid

:class:`MomentumGrid` object, containing a grid of momentum points

.. attribute:: timeGrid

:class:`TimeGrid` object containing timestep vector amongst other


Parameters for autoInitialGrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: NxiScalingFactor

uniformly rescales the predicted Nxi value from 
autoInitialGrid by a factor of NxiScalingFactor (default 1)

.. attribute:: dyBulkScalingFactor

uniformly rescales the desired grid spacing dy
at y = 0 used by autoInitialGrid (default 1, lower value yields higher resolution)

.. attribute:: dyTailScalingFactor

uniformly rescales the desired grid spacing dy
at y = yMax used by autoInitialGrid (default 1, lower value yields higher resolution)

.. attribute:: Nxi_min

.. attribute:: Nxi_max

.. attribute:: pMax_ceiling

.. attribute:: pMaxIncreaseFactor

.. attribute:: minPMaxMarginFactor

.. attribute:: pSwitch

.. attribute:: percentBulk

.. attribute:: r

Functions
--------------

Constructor
%%%%%%%%%%%%%%%%%

.. function:: this = State(TRef,nRef)

Construct a new state with :class:`MomentumGrid`, :class:`TimeGrid`, :class:`PhysicalParams` and :class:`Reference` classes.

.. function:: setInitialRunaway(this, Distribution)

Sets an initial runaway current to be used for autinitial grid.
Distribution is :class:`Distribution` object from which the current is calculated from.
The autInitGrid then takes into account for already created runaways when estimating yMax and other parameters.

.. function:: autoInitGrid(this,useScreening,useInelastic)

Automatically sets yMax, Nxi, Ny and gridWidth for gridMode 6 given the physical scenario. 
Unless Initial runaway is set, theory using a gaussian distribution as start distribution is used to estimate how far the tail will reach.
Requiers that all other physical and time parameters are set to their values to give meaningfull results.
