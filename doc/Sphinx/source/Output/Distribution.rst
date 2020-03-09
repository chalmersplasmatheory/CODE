Distributions
=================

.. classdef:: Distribution 

Properties
---------------

Distribution function
%%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: f 

Disitrubition function f 

.. math::
        f = \begin{bmatrix}
            \begin{bmatrix} 
                f_1(y)\\
            \end{bmatrix}\\
            \begin{bmatrix} 
                f_2(y)\\
            \end{bmatrix}\\
            \vdots\\
            \begin{bmatrix} 
                f_{N_\xi}(y)\\
            \end{bmatrix}
        \end{bmatrix}

where f_i(y) is the projection on the i_th Legandre mode at momentum y where y is defined in :attr:`momentumGrid`.

.. attribute:: momentumGrid

:class:`MomentumGrid` where the distribution is defined

Design note, momentum grid is not checked for changes when saved to
this class. Therefore if things change in memory to which this points
to, it will be unnoticed and wrong.

.. attribute:: time

Time where distribution is defined 

Reference Values
%%%%%%%%%%%%%%%%%%%%

reference values as :class:`Reference` object in :attr:`momentumGrid`. 
Just saved another time. 
This will prob dissapear in future version. 
Access with [Disitribution].momentumGrid.[TRef,nRef,...] instead

.. attribute:: TRef

.. attribute:: nRef

.. attribute:: deltaRef

.. attribute:: lnLambdaRef

.. attribute:: nueeRef


Functions
----------------

Constructor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. function::  Distribution(f,momentumGrid,time)

Constructs a new Distribution
