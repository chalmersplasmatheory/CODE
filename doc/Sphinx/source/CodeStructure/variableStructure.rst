Code description
======================

In this chapter we describe the overall structure of how indexing is done in various variables

Geneneral time concept
----------------------

In variables in one dimension that depends on time, the indexing start on 1 which corresponds to the values at the initial distribution (whether the initial distribution is given or not).
In variables with more than one dimension, one must check the specific variables documentation to know in which dimension (row/col) the time is defined.
Sometimes there are also variables which have two dimensions which of none represent time advancment.
Such can for example be when plotting the distribution function i 2D momentum space.

Structure of Common Variables
------------------------------
Here the some variables structure is define.

The distribution function is a 1d array (and is for example present in :class:`SmartLUSolver` as variable fatTimeIndex or in :class:`Distribution` as variable :attr:`f`). 
Each coloumn has the structure

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

where f_i(y) is the projection on the i_th Legandre mode at momentum y.
The associated momentum vector is found in the :class:`State` or :class:`MomentumGrid` objects which are properties to the class.

Time dependent variables. 
In physical params class, :attr:`E`, :attr:`T`, :attr:`n` and :attr:`Z` can be time dependent.
From the users perspective, these are always defined at the associated :class:`TimeGrid` (unless half complicated manouvers are done, see later in this documentation).
The raw input data is always saved and it is from this data the interpolation is done.

Equation 
-------------------------

The equation beeing solved is the kinetic equation, which reads 

.. math::

        \frac{df}{dt} + \vec{F} \cdot \frac{df}{d\vec{p}} + \vec{v} \cdot \frac{df}{d\vec{x}} = C \{ f \}.

In this code, only momentum is considered and all spatial dependence is disregarded. 
The momentum is parameterized using spherical coordinates, where the polar angle is replaced with the pitch angle.
The distribution function is assumed to be independent of the azimuthal angle due to gyro averaging.
The pitch angle is described by its cosine as

.. math::

        \xi = \cos (\theta)

where theta is the pitch angle.

By only considering the average electric field along the magnetic field the equation

.. math::

        \frac{df}{dt} + e E (\xi \frac{d f}{d p} + \frac{1 - \xi^2}{p} \frac{d f}{d \xi})
        = C\{f\}

is obtained; where p is the momentum, E the electric field strength and e the electron charge. 
This is then discritized as 

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

where f_i(y) is the projection on the i_th Legandre mode at momentum y.

Numerically this is then rewritten as a matrix equation

.. math:: 

       \overleftrightarrow{O}_{\text{Implicit}} f^{i+1} + \frac{f^{i+1}-f^{i}}{\Delta t} = \overleftrightarrow{O}_{\text{Explicit}} f^{i} + S(f).

Operators extending :class:`ImplicitOperator` is written to be a part of the operator indexed ''Implicit'' and operators extending :class:`ExplicitOperator` is written to be part of the operator indexed Explicit.
The classes extending the :class:`Source` is written to be added to the S(f) function and is done so since their matrix representation is large. 
Therefore it is more numerically feasible to handle them explicitly as functions of f.
The difference between :class:`Source` and :class:`ExplicitOperator` is that :class:`Source` uses the distribution function all the time and recalculates the source whilst the :class:`ExplicitOperator` calculates a matrix representation and then uses the matrix product to calculate the term.

The equation is normalized to a reference temperature and a reference density such that the unit of time is in electron electron collision times using the most probable thermal speed in the collision frequency. 
The distribution function is normalized such that a 1D-maxwellian distribution with this normalization and density of reference density acquire a value of 1 at velocity zero.
The full derivations and equations is described in ''Reference To Adams Paper''. 
