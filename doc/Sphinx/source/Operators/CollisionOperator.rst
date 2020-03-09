CollisionOperator
=================
.. class:: CollisionOperator

extends :class:`ImplicitOperator`.
Class controlling the collsion operator Operator of the equation system. Right now gives support to 4 collision operators. One fully linearized, one fully relativistic (with relativistic bulk) and one with relativistic fast particles but with classical maxwellian bulk. ``insert info about all the collision operators``.

Functions
------------

.. function:: generateOperatorMatrix(this,iteration)

Generates the Operator matrix for the iteration, physical quantities taken from the plasma class. Matrix is valid for the grid defined in grid class (number of Legendre Modes, number of radial grid points and exact spacing of the grid points).

.. function:: CalculateEnergyDependentlnLambda(this, p, y, lnLambda0)

returns quantities relevant for Energy dependent lambda. 
This is an internal function and therefore returns are normalized to fit in with rest of normalization of the code.
Calculates the enhancement of the Coulomb logarithm replacing the termal speed with the relativistic momentum at high energy. 
Basically, nuS -> nuS*(lnLambdaHat + Ghat). 
The formula is a simple interpolation between the energy dependent expression and the thermal speed expression at low energy. 
References:
Wesson p. 792 and (lnL0 = lnL^ee) and Solodev-Betti (2008) for energy dependence

OUTPUT: 

* lnLambdaeehat -
* dlnLambdaeehatdp - derivative of lnLambda with respect to momentum
* lnLambdaeihat -
* boltzCorrectionlnLambdaeeHat - 
* boltzCorrectiondlnLambdaeeHatdp - 

All have the same size as momentum grid.

.. function:: CalculateInelasticEnhancement(this, species, p, y,lnLambda0, useInelastic, iteration )

calculates the enhancement of the
electron-electron slowing down frequency (\nu_S^{ee}) due to inelastic
collisions with bound electrons around the ion. Basically, nuS->
nuS\*(lnLambdaHat + Ghat). 

INPUTS:

* species: struct with fields:
        * nj(jspecies, times) - matrix containing number density of each species as a function of time (m^{-3})
        * ZAtomicNumber - atomic number of each species. row vector, not a function of time
        * Z0NetCharge - net charge for each species. row vector, not a function of time
        * times - time steps. If constant, just set i to 1.

* p = normalized momentum gamma*m*v/(m*c) = delta*y
* LnLambda the coulomb logarithm
* useInelastic - which model to use
        0: Ignore inelastic collisions
        1: Bethe formula with values of mean excitation energy from Sauer et al.
        2: RP Rule of thumb (cutoff at y=10 to conserve particles)


OUTPUT: 

* Ghat - 
* dGhatdp - derivative of Ghat

Both have the same size as y.

.. function:: CalculateSpeciesScreeningFactor(this, useScreening,Z,Ne,Z0,p,Lmax )

CalculateSpeciesScreeningFactor calculates the function g(p): the \nu =
(1+g/lnLambda), for each legendre mode. p must be a row vector.
output: g, with size (Lmax, Ny). 

INPUTS:

* species - with the fields
        * nj(jspecies, times) - matrix containing number density of each
        * species as a function of time (m^{-3})
        * ZAtomicNumber atomic number of each species. row vector, not a function of time
        * Z0NetCharge net charge for each species. row vector, not a function of time
        * times time steps. If constant, just set i to 1.

* p = normalized momentum gamma*m*v/(m*c) = delta*y
* LnLambda the coulomb logarithm. It is now a (Nxi-1, Ny ) matrix! (p-dependent)
* useScreening: which model to use
        0. Ignore screening (complete screening)
        1. Fokker--Planck collision operator, TF model (works for all species)
        2. Fokker--Planck collision operator, TF-DFT model (works for some ions)
        3. Fokker--Planck collision operator, full DFT model (works for some ions)
        4. Boltzmann collision operator, TF model
        5. Boltzmann collision operator, TF-DFT model
        6. Boltzmann collision operator, full DFT model
        7. no screening = full penetration

OUTPUT: 

* g - 

.. function:: CalculateScreeningFactor(this, species, p, Lmax, LnLambda, useScreening, iteration)

CALCULATESCREENINGFACTOR calculates the enhancement of the electron-ion
deflection frequency from screening effects

INPUTS:

* species - struct with the fields
        * nj(jspecies, times) - matrix containing number density of each
        * species as a function of time (m^{-3})
        * ZAtomicNumber - atomic number of each species. row vector, not a function of time
        * Z0NetCharge - net charge for each species. row vector, not a function of time
        * times - time steps. If constant, just set i to 1.

*  p = normalized momentum gamma*m*v/(m*c) = delta*y
*  LnLambda the coulomb logarithm. It is now a (Nxi-1, Ny ) matrix! (p-dependent)
*  useScreening: which model to use
	0. Ignore screening (complete screening)
	1. Fokker--Planck collision operator, TF model (works for all species)
	2. Fokker--Planck collision operator, TF-DFT model (works for some ions)
	3. Fokker--Planck collision operator, full DFT model (works for some ions)
	4. Boltzmann collision operator, TF model
	5. Boltzmann collision operator, TF-DFT model
	6. Boltzmann collision operator, full DFT model
	7. no screening = full penetration


OUTPUT: 

* screeningFactor - is a matrix of size (Nxi, Ny), that describes the enhancement of the deflection frequcency for the each Legendre mode and all y values on the grid.

.. function:: CalcgBoltzmann(~, F,Lmax, Z0, Ne )

calculates the Legendre-mode-dependent correction to the
Fokker-Planck deflection frequency neglecting screening.

INPUT:

* F(p) -  a function handle (if numeric, pchip(F,p) works
* k -  p/mc (vector)
* Z0 - net charge
* Ne: number of bound electrons
* Lmax - maximum legendre mode

The integrals are split at \theta_cutoff because the legendres behave badly for too small arguments

OUTPUT:

* boltzmanng - 
* p -
