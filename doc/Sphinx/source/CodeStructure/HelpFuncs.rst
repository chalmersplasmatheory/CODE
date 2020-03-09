cs = CalculateMollerCS(g,g1)
-------------------------------------------------------------------

.. function:: cs = CalculateMollerCS(g,g1)

Calculates the Møller cross-section (cs) for large-angle
electron-electron collisions between an incoming electron with initial
:math:`gamma=g1` (relatistic gamma function) and an electron at rest, which
acquires :math:`gamma = g` in the collision. 
Input g and g1 can be arrays (of the same size).

moment = CalculateRunawayMomentXiDependent(Nxi,Ny,f,xiIntMat,weightsMat)
------------------------------------------------------------------------------

.. function:: moment = CalculateRunawayMomentXiDependent(Nxi,Ny,f,xiIntMat,weightsMat)

Calculates the :math:`\xi` dependent runaway moment given integration matricices.

Input:

    Nxi - number of legandre modes

    Ny - Number of grid points

    f - distribution function

    xiIntMat - 

    weightsMat -

sigma = CalculateSigmaIntegral(gamma,gamma_m)
-------------------------------------------------------------------

.. function:: sigma = CalculateSigmaIntegral(gamma,gamma_m)

Calculates the integral

.. math::

    \int_{\gamma_m}^{\gamma+1-\gamma_m} \Sigma(\gamma,\gamma_1) d\gamma_1

analytically. 
This calculation is needed for the full, conservative close-collision operator.
Gamma can be a vector or array, :math:`\gamma_m` must be a scalar. 
If :math:`\gamma+1-\gamma_m< \gamma_m`, the integral is set to be zero.


[m,xiGrid] = GenerateCumIntMat(Nxi,nPoints)
------------------------

.. function:: [m,xiGrid] = GenerateCumIntMat(Nxi,nPoints)

Generates grids to cumulatively integrate legendre modes over a xi grid

f = interpolateDistribution(dist,MG)
-------------------------------------------------------------------

.. function:: f = interpolateDistribution(dist,MG)

Interpolates each legandre mode function from distribution for itself with the momentumGrid in query.
The distribution in Distribution is interpreted to be zero above its maximum momentum, and in legandre modes higher than its Nxi.
Input:

    dist - :class:`Distribution` object, from which the distribution should be interpolated.
    
    MG - :class:`MomentumGrid` object, to which grid the output distribution is interpolated.

outPls = LegendrePolynomials(l,x)
-------------------------------------------------------------------

.. function:: outPls = LegendrePolynomials(l,x)

LegendrePolynomials calculates :math:`P_l(x)`.

Calculates the legendre polynomials :math:`P_i(x)` for :math:`i=0,1,...,l` using
Bonnet's recursion formula:

.. math::

    (n+1)*P_{n+1}(x) = (2n+1)*x*P_n(x) - n*P_{n-1}(x),

where :math:`P_0(x) = 1` and :math:`P_1(x) = x`.

Usage:
pls = LegendrePolynomials(l,x)

l is the (highest) mode number and x is the coordinate, which must be
a row vector (not a matrix). 

pls has the structure

.. math::

    \begin{matrix}
    [ P_0(x) ]\\
    [ P_1(x) ]\\
    [   ...  ]\\
    [ P_l(x) ]
    \end{matrix}


[x,w] = lgwt(N,a,b)
-------------------------------------------------------------------

.. function:: [x,w] = lgwt(N,a,b)

This script is for computing definite integrals using Legendre-Gauss
Quadrature. Computes the Legendre-Gauss nodes and weights on an interval
[a,b] with truncation order N

Suppose you have a continuous function f(x) which is defined on [a,b]
which you can evaluate at any x in [a,b]. Simply evaluate it at all of
the values contained in the x vector to obtain a vector f. Then compute
the definite integral using sum(f.*w).

Written by Greg von Winckel - 02/25/2004

[x, w, D, DD] = m20121125_04_DifferentiationMatricesForUniformGrid(N, xMin, xMax, scheme)
-------------------------------------------------------------------------------------------------

.. function:: [x, w, D, DD] = m20121125_04_DifferentiationMatricesForUniformGrid(N, xMin, xMax, scheme)

Finite difference differentiation matrices and integration weights for a
uniform grid.

Created by Matt Landreman,
Massachusetts Institute of Technology, Plasma Science & Fusion Center, 2012.

Inputs:

    N - number of grid points.

    xMin - minimum value in the domain.

    xMax - maximum value in the domain.

    scheme - switch for controlling order of accuracy for differentiation and handling of endpoints.

Options for scheme:

    0 -  The domain [xMin, xMax] is assumed to be periodic. A 3-point stencil
    is used everywhere. A grid point will be placed at xMin but not
    xMax.

    1 -  Same as scheme=0, except a grid point will be placed at xMax but not
    xMin.

    2 -  The domain [xMin, xMax] is assumed to be non-periodic. A 3-point
    stencil is used everywhere.  The first and last row of the
    differentiation matrices will use one-sided differences, so they
    will each have a non-tridiagonal element.

    3 -  The same as scheme=2, except that the first differentiation matrix
    will use a 2-point 1-sided stencil for the first and last elements
    so the matrix is strictly tri-diagonal.  The 2nd derivative matrix
    is the same as for option 2, since it is not possible to compute
    the 2nd derivative with only a 2-point stencil.

    10 - The domain [xMin, xMax] is assumed to be periodic. A 5-point stencil
    is used everywhere. A grid point will be placed at xMin but not
    xMax.  This option is like scheme=0 but more accurate.

    11 - Same as scheme=10, except a grid point will be placed at xMax but
    not xMin.  This option is like scheme=1 but more accurate.

    12 - The domain [xMin, xMax] is assumed to be non-periodic. A 5-point
    stencil is used everywhere.  The first two and last two rows of
    the differentiation matrices will then each have non-pentadiagonal
    elements.

    13 - The same as option 12, except that 3-point stencils are used for the
    first and last rows of the differentiation matrices, and 4-point
    stencils are used for the 2nd and penultimate rows of the
    differentiation matrices.  With this option, both differentiation
    matrices are strictly penta-diagonal.

    20 - The domain [xMin, xMax] is assumed to be periodic. Spectral
    differentiation matrices are returned. A grid point will be placed
    at xMin but not xMax.

    21 - Same as scheme=20, except a grid point will be placed at xMax but not
    xMin.

Outputs:

    x - column vector with the grid points.

    w - column vector with the weights for integration using the trapezoid rule.

    D - matrix for differentiation.

    DD - matrix for the 2nd derivative.

m = MapToXiIntMat(cumintMat,xiGrid,EOverEc,y2,deltaRef)
------------------------------------------------------------------------

.. function:: m = MapToXiIntMat(cumintMat,xiGrid,EOverEc,y2,deltaRef)

TODO Documenation


[posF,negF] = SumLegModesAtXi1(f,Ny,Nxi)
-------------------------------------------------------------------

.. function:: [posF,negF] = SumLegModesAtXi1(f,Ny,Nxi)


Calculates the total distribution at xi=1 (v_perp=0) by summing
over all the Legendre modes of the distribution

Input:

    f - distribution function with normalization as in CODE
    and Nxi legandre modes projections. first Ny points are for first
    legandre mode.

    Nxi - Number of Legandre Modes

    Ny - Number of radial grid points


i = CalculateCurrentAlongB(f,y,yWeights, Ny, nRef ,deltaRef)
-------------------------------------------------------------------

.. function:: i = CalculateCurrentAlongB(f,y,yWeights, Ny, nRef ,deltaRef)


Calculates current density in SI along magnetic axis
specifically calculates 

.. math::

    n \iiint \hat{f} e v d^3v

which reduces to

.. math::

    4 n_{Ref} e v_{Ref} / (3 \sqrt(\pi) m_e^3) \int_0^\infty y^3 / \gamma dy


where :math:`e` is electron charge, and :math:`v` is the electron speed

Input:

    f - distribution %function with normalization as in CODE
    and Nxi legandre modes projections. first Ny points are for first
    legandre mode.

    y - is vector of momentum points as y = gamma v / v_Ref where gamma is
    relativistic gamma function, v electron speed, and v_Ref reference
    thermal speed of electrons used in y

    yWeights - contains dy weights for integrating over y

    nRef - referece density of electrons in SI

    deltaRef - vRef/c where c is speed of light

n = CalculateDensity(f,Ny,Nxi,y,yWeights,nRef)
-------------------------------------------------------------------

.. function:: n = CalculateDensity(f,Ny,Nxi,y,yWeights,nRef)

Calculates density n as integral over momentum space. The density is output in SI.

Input:

    f - distribution %function normalized as in CODE

    y - momentum grid

    yWeights - weights for integrating y

    nRef - reference density

    Nxi - Legendre Modes

    Ny - points in y;

delta = getDeltaFromT(T)
-------------------------------------------------------------------

.. function:: delta = getDeltaFromT(T)

Computes the delta (=v_th/c) parameter used in CODE from the
electron temperature (in eV).

[delta,lnLambda,collFreq,BHat] = getDerivedParameters(T,density,B)
-------------------------------------------------------------------

.. function:: [delta,lnLambda,collFreq,BHat] = getDerivedParameters(T,density,B)

Calculates the parameters delta, coloumb logarithm lnLambda and the
collision frequency (and time) for given E (in V/m), T (in eV) and density
(in m^{-3}). Collission frequency uses Matt's definition
Collision time = 1/Collision Frequency

[EOverEC,EOverED,EHat] = getNormalizedEFields(E,T,density)
-------------------------------------------------------------------

.. function:: [EOverEC,EOverED,EHat] = getNormalizedEFields(E,T,density)

Normalizes an electric field in V/m to the critical field EOverEc,
and to the dreicer field EOverED
E is the field in V/m, T is the electron temperature in eV and
density is the electron density in m^{-3}.

sigma = GetSpitzerConductivity(T,n,Z)
-------------------------------------------------------------------

.. function:: sigma = GetSpitzerConductivity(T,n,Z)

Calculates the plasma Spitzer conductivity in 
units of (Ohm m)^-1, from the formula on p. 72 of Helander & Sigmar. 
Usage:
    
sigma = GetSpitzerConductivity(T,n,Z)

n must be in units of m^(-3), T in eV.
There is a prefactor that depends on Z which is only known numerically
(table on p.74 in Helander & Sigmar and Table III in Spitzer & Härm).
Let's interpolate to find the value for any Z. Since there is a data
point at infinity, lets use 1/Z for the interpolation.
Simplified expression in Chang (Eq. [5-76]):
sigma = 19.2307e3*T.^(1.5)./Z./lnLambda;
Written by Adam Stahl
