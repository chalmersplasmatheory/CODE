TimeGrid 
============

Class describing and containing the descritization in time.

.. class:: TimeGrid

Properties
-----------------

Time advance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: dt

timestep size (normalized). Distance between different time steps in uniform grid. In case of nonuniform, determines the order of magnitude of timedifference

.. attribute:: tMax

maximum time (minimum is so far always 0) (normalized). Maximum time to simulate to. This will always be exact (to machine precision)

.. attribute:: timeStepMode

specifies how the time step vector is generated

0. Constant time step (dt), in units of timeUnit

2. Logarithmic -- use progressively longer time steps. Usefull
for convergence towards a steady state. dt is used for the
first time step.

3. Stepwise logarithmic -- like 2, except that several time
steps are taken with each step length. This is to avoid
rebuilding the matrix in every time step.

.. attribute:: logGridScaling

step size scaling for the logarithmic grid

.. attribute:: logGridSubSteps

how many steps to take for each time step length
with timeStepMode 3

.. attribute:: logGridMaxStep

maximum time step allowed in timeStepModes 2 and 3

.. attribute:: timesteps

actual times the distribution is calculated in (normalized)

.. attribute:: dts

vector of all small timechanges in (normalized)

.. attribute:: dtsHasChanged

vector contantaining if the a element in dts vector
has changed or not

.. attribute:: nTimeSteps

number of times where the distribution is defined
(actually one greater than the number of timesteps
that should be taken)

Misc
%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: PhysicalParams

:class:`PhysicalParams` object.

Design note: intention is that physicalParameters containing an instance of
this class also is contained in this object, think of it as pairs. 
Right now it is possible to 'hack' the construction by creating a
TimeGrid object, and a PhysicalParams. Then setting the TimeGrid in
the PhysicalParams. Afterwards it is possible to create a new 
PhysicalParams and set the already
created TimeGrid in the new PhysicalParams. The firstly created PhysicalParams 
now has a TimeGrid object in it with a PhysicalParams object in it which 
is not pointing to itself. The 'pair' structure is then broken. Two
fixes are availible: one seperating the TimeGrid to a copy of itself
(new pointer) and using the copy for the old TimeGrid or save both
TimeGrid objects in the same reference object.
properties (SetAccess = protected)
physicalParams 

.. attribute:: reference

:class:`Reference` object shared amongst all state objects, containing reference values.

Functions
----------------

Constructor
%%%%%%%%%%%%

.. function:: TimeGrid(reference,varargin)


Set functions
%%%%%%%%%%%%%%%%%%%

.. function:: setPhysicalParams(this,phP)

Sets phP to passed value and also reinitializes the time grid

.. function:: setResolution(this,varargin)

Sets varargin to passed value and also reinitializes the time grid

.. function:: setdt(this, dt)

Sets dt to passed value and also reinitializes the time grid

.. function:: settMax(this, tMax)

Sets tMax to passed value and also reinitializes the time grid

.. function:: settimeStepMode(this, timeStepMode)

Sets timeStepMode to passed value and also reinitializes the time grid

.. function:: setlogGridScaling(this, logGridScaling)

Sets logGridScaling to passed value and also reinitializes the time grid

.. function:: setlogGridSubSteps(this, logGridSubSteps)

Sets logGridSubSteps to passed value and also reinitializes the time grid

.. function:: setlogGridMaxStep(this, logGridMaxStep)

Sets logGridMaxStep to passed value and also reinitializes the time grid

.. function:: initializeTimeGrid(this)

.. function:: getdt(this,timeUnit)

Returns dt in specified time unit

Input:

    timeUnit
        
        s - seconds
    
        ms - milliseconds

        normalized - in nueeRef from 


.. function:: gettMax(this,timeUnit)

Returns tMax in specified time unit

Input:

    timeUnit
        
        s - seconds
    
        ms - milliseconds

        normalized - in nueeRef from 


.. function:: gettimeStepMode(this)

Returns timeStepMode 

Returns timeStepMode 

.. function:: getlogGridScaling(this)

Returns logGridScaling 

Returns logGridScaling 

.. function:: getlogGridSubSteps(this)

Returns logGridSubSteps 

Returns logGridSubSteps 

.. function:: getlogGridMaxStep(this)

Returns logGridMaxStep 

Returns logGridMaxStep 

.. function:: gettimesteps(this,timeUnit)

Returns timesteps in specified time unit

Input:

    timeUnit
        
        s - seconds
    
        ms - milliseconds

        normalized - in nueeRef from 

.. function:: getdts(this,pts,timeUnit)

Returns dts in specified timeUnit

Input:

    timeUnit
        
        s - seconds
    
        ms - milliseconds

        normalized - in nueeRef from 

.. function:: getdtsHasChanged(this)

Returns dtsHasChanged 

.. function:: getnTimeSteps(this)

Returns nTimeSteps 

