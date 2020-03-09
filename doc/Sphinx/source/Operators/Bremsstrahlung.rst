Bremsstrahlung
===================

.. class:: Bremsstrahlung < ExplicitOperator

BOLTZMANN generates the bremstrahlung matrix 
basically does the things GenerateMatricesForBoltzmannOperators did
minus the knock on matrix plus doing the line:
boltzmannMatrix = nBar*(1+this.plasma.Z(iteration))*bremsMatrix

Properties
-----------------------

.. attribute:: NxiUsed

What number of legandre modes used.

.. attribute:: pUsed

What momentum vector used

.. attribute:: NyInterpUsed

What number interpolation points of y used.

.. attribute:: useScreeningUsed

What screening switch used in :class:`EquationSettings`.

.. attribute:: bremsModeUsed

What bremsmode is calculated

.. attribute:: speciesUsed

Used species. Todo: smart implementation to check that only the current timestep is checked if same

.. attribute:: nBarUsed

What density used

.. attribute:: ZUsed

Charge used

.. attribute:: deltaUsed

Velocity over speed of light used.

.. attribute:: lnLambdaRefUsed

Reference coloumb logarithm used.

Functions
----------------
    
.. function:: this = Bremsstrahlung(state,eqSettings)

Construct a new instance of this class

.. function:: matrixHasChanged = generateOperatorMatrix(this,iteration)

Call this function to generate the matrix.

.. function:: brMat = GenerateBremsMatrix(this,iteration)

.. function:: [M,M_source,M_sink] = GenerateBremsstrahlungMatrix(this,mode,smallK)

Helpfunctions for building the matricies 
************************************************

.. function:: d1sigma= d1sigma_BetheHeitler(~,p,k)

.. function:: d1sigma_screened = d1sigma_BetheHeitler_screened(this,p0,k,species)

.. function:: C = GeneratePitchOperator(~,p,Nxi,kMin,kMax)

.. function:: [ ZMinusFSquared ] = calculateAveragedFormFactorBrems(this, species, iteration)

.. function:: W = diffCrossSec(this,p,p1,theta)

.. function:: outPl = LastLegendrePolynomial(this,l,x)

.. function:: K = BuildInterpolationMatrix(this,x,xp)

