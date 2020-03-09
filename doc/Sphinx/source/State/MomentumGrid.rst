MomentumGrid
======================
.. class:: MomentumGrid

Class for the grid in momentum space.

Properties
---------------------

Reference values
%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: reference

:class:`Reference` object

Resolution parameters:
%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: Nxi

number of Legendre modes used to represent the pitch-angle
coordinate (determines resolution in xi=cos(theta))

.. attribute:: Ny

number of points in the grid used to represent momentum

.. attribute:: yMax

maximum of the momentum grid, in units of y=gamma*v/v_{th}.
Together, Ny and yMax determine the resolution in momentum

.. attribute:: y

momentum grid p/p_thermal, where p_thermal is treated
classically

.. attribute:: y2

y.^2

.. attribute:: y4

y.^4

.. attribute:: x

v/v_thermal
x2. x.^2

.. attribute:: gamma

Relativistic gamma function

Differential Calculation Params: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: ddy

Differential operator, differenting with respect to y, requiers
vector of length same as y

.. attribute:: d2dy2 

Same as ddy but second derivetive

.. attribute:: yWeights

Integration operator with respect to y, same limitation as
ddy (row vector)

.. attribute:: ddx

Differential operator with respect to x
d2dx2 - second order differential operator with respect to x

Numerical grid parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: yGridMode

Specifies how the points on the momentum grid are chosen

0. uniform grid on [0, dy*(Ny-1)]

1. nonuniform grid, remapped using y = scaleFactor*tan(pi*a*s/2)

2. nonuniform grid, remapped using y = -scaleFactor*ln(a-s)

3. nonuniform grid, remapped using y = scaleFactor*s/(a-s)

4. nonuniform grid, remapped using y = s^2 + gridParameter*s,
where gridParameter>0 (default)

5. nonuniform grid, smooth tanh step between dense bulk and

sparse tail. gridParameter controls the spacing in
the bulk (and is in units of y). gridStepWidth
controls the width of the smooth step between bulk
and tail, and gridStepPosition its position (in
units of y_c). gridParameter=1/15,
gridStepWidth=1/50 and gridStepPosition=2 is a good
starting point, not implemented yet
Above, s is a uniform grid and a is a constant close to 1.
When running a CODE+GO calculation, yGridMode must be
either 0 or 4. The reason is that only these two remappings
are able to handle varying grids in the required way.

.. attribute:: gridParameter

for use with yGridMode 4 and 5

.. attribute:: gridStepWidth

for use with yGridMode 5

.. attribute:: gridStepPosition

for use with yGridMode 5

.. attribute:: yMaxBoundaryCondition

Boundary condition at yMax (does this more belong in the equation settings?)

1. Dirichlet: F=0 (default)

2. Robin: dF/dy + (2/y)*F=0, which forces F to behave like 1/y^2

3. Do not apply a boundary condition at yMax

4. Dirichlet: F=0, with a bit of extra d2dy2 added at the last few grid points to eliminate ringing

.. attribute:: artificialDissipationStrength

only used with yMaxBoundaryCondition=4
to control the amount of ringing

.. attribute:: artificialDissipationWidth

in momentum (how close to yMax). Only
used with yMaxBoundaryCondition=4 to
control the amount of ringing


.. attribute:: VERSION 

Functions
-------------------

Constructor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. function:: this = MomentumGrid(reference,varargin)

creates new momentum grid

Set functions
%%%%%%%%%%%%%%%%%%%%%

.. function:: setResolution(this,varargin)

Sets resolution parameters, standard matlab syntax ('name',value). Also reinitializes the momentum grid

.. function:: setNxi(this, Nxi)

Sets Nxi to passed value and also reinitializes the momentum grid

.. function:: setNy(this, Ny)

Sets Ny to passed value and also reinitializes the momentum grid

.. function:: setyMax(this, yMax)

Sets yMax to passed value and also reinitializes the momentum grid

.. function:: setyGridMode(this, yGridMode)

Sets yGridMode to passed value and also reinitializes the momentum grid

.. function:: setgridParameter(this, gridParameter)

Sets gridParameter to passed value and also reinitializes the momentum grid

.. function:: setgridStepWidth(this, gridStepWidth)

Sets gridStepWidth to passed value and also reinitializes the momentum grid

.. function:: setgridStepPosition(this, gridStepPosition)

Sets gridStepPosition to passed value and also reinitializes the momentum grid

.. function:: setyMaxBoundaryCondition(this, yMaxBoundaryCondition)

Sets yMaxBoundaryCondition to passed value and also reinitializes the momentum grid

.. function:: setartificialDissipationStrength(this, artificialDissipationStrength)

Sets artificialDissipationStrength to passed value and also reinitializes the momentum grid

.. function:: setartificialDissipationWidth(this, artificialDissipationWidth)

Sets artificialDissipationWidth to passed value and also reinitializes the momentum grid

Misc.
%%%%%%%%%%%%%

.. function:: initializeMomentumGrid(this)

Creates momentum vectors, derivatives and more
