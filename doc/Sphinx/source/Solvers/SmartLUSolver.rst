SmartLUSolver
=========================

.. class:: SmartLUSolver

Extends :class:`Solver` class

Properties
------------------

.. attribute:: timeIndex

Time index where simulation will start and where :attr:`fattimeIndex` is defined.

.. attribute:: fattimeIndex

Distribution function that will be used to simulate next timestep and is thus defined at :attr:`timeIndex`.

.. attribute:: nSaves

How many saves to do when invoking :func:`takeTimeSteps`

.. attribute:: saveIndices

Which indecies will be saved.
Initialized after nSaves have been set and is set when invoking :func:`takeTimeSteps` with :func:`initsaveIndices`. 

.. attribute:: factor_L

L factor in LU factorization for last taken step in time.

.. attribute:: factor_U

U factor in LU factorization for last taken step in time.

.. attribute:: factor_P

P factor in LU factorization for last taken step in time.

.. attribute:: factor_Q

Q factor in LU factorization for last taken step in time.

Functions 
--------------------

.. function:: this = SmartLUSolver(state,eqSettings)

Construct a new instance of this class
    
.. function::  output =  takeTimeSteps(this,output)
        
Takes the remaining timesteps, starts from :attr:`timeIndex` and ends at end of :class:`TimeGrid` object in :class:`State` object in this class.

If output is supplied, this will be modified by adding new output to it, usefull when
extending the distributions.
Remember that this is pointers, so you if you supply a output class you only need to write takeTimeSteps(output); (e.i. no assignment as it is a handle class).

Setters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. function::  setf(this,distribution,momentumGrid)
        
Sets the fattimeIndex to the distribution (is a :class:`Distribution` object).
If momentumGrid (:class:`MomentumGrid` object) is supplied then it interpolates the distribution function from distribution to the supplied momentumGrid object.

.. function::  setTime(this,time,timeUnit)
       
Sets the :attr:`timeIndex` to last index where time is greater than the :attr:`timesteps` in :class:`TimeGrid` object in :attr:`state`.

.. function::  setnSaves(this,nSaves)
   
Set how many steps will be returned.
Must be greater than 2.

initializing functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. function::  initsaveIndices(this)

Sets what indices to save. 
Is called by :class:`takeTimeSteps`.
Always saves first and last index
        
.. function::  clearf(this)
        
Reset the distribution function :attr:`fattimeIndex` to be empty.

.. function::  setftoZero(this)
        
Set the distribution function :attr:`fattimeIndex` to all zeroes.

.. function::  resettimeIndex(this)
        
Sets the :attr:`timeIndex` to be one, e.i. restart the simulation from first timestep.

.. function::  setftoMaxmellian(this)

Sets the distribution function :attr:`fattimeIndex` to be Maxmellian.
        
.. function::  rhs0 = enforceParticleAndHeat(this)
            
this adds a particle and heat source dependent on what settings are
used in eqSetings. These are added to the 0th Legandre Mode.

