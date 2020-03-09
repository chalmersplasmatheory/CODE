EquationSettings
==================
.. class:: EquationSettings

.. attribute:: operators

Operators, in a cell array, necassary to fullfill the settings specified below.

.. attribute:: useScreening

take screening into account

0. Ignore screening (complete screening, correct at low energy)
1. Fokker--Planck collision operator, TF model (works for all species)
2. Fokker--Planck collision operator, TF-DFT model (works for some ionization degrees of argon and Xe)
3. Fokker--Planck collision operator, full DFT model (only works for Ar+ so far)
4. Boltzmann collision operator, TF model
5. Boltzmann collision operator, TF-DFT model
6. Boltzmann collision operator, full DFT model
7. Full penetration (no screening, Z0^2-> Z^2, correct in the limit p-> infinity)

.. attribute:: useInelastictake 

inelastic collisions with bound electrons intoaccount 

0. Ignore inelastic collisions 
1. Bethe formula 
2. RP

.. attribute:: useEnergyDependentLnLambdaScreening

Logarithmic enhancement of the Coulomb logarithm 

0. No (standard formula)
1. Yes (change in the collision operator). Interpolated Landau form.
2. Enhance by an ad-hoc factor of 1.3.

Collision operators and conservational properties
---------------------------------------------------


.. attribute:: collisionOperator 

specifies the type of collision operator to use 

0. Fully relativistic test-particle collision operator. NOT momentum-conserving. Relativistic Fokker-Planck coefficients from OJ Pike & SJ Rose, Phys Rev E 89 (2014). 
1. Non-relativistic, linearized, momentum-conserving. Includes modified Rosenbluth potentials. 
2. Use the non-relativistic field particle term from 1, togetherwith the relativistic test particle term from 4. Generally works better than 1, as the tail quickly becomes relativistic. The field-particle term only affects the  bulk, so the non-relativistic approximation is usually good  there.
3. Use an ad-hoc relativistic generalization of the field-particle term from 1, together with the test-particle terms from 4. Generally appears to work better than 1 & 2,  although the physics basis is shaky. !!!Experimental!!! 
4. Relativistic collision operator for non-relativisticthermal speeds. Agrees with 0 when T << 10 keV. See the CODE-paper for details.

.. attribute:: useNonRelativisticCollOp 

usethe non-relativistic limit of thecollision operator, for instance forbenchmarking hot-tail generation againstSmith et al. (2005), or Smith andVerwichte (2008).

Avalanche (secondary runaway) source term
-------------------------------------------

.. attribute:: sourceMode 

the type of avalanche source to use 

0. do not include a source
1. include a Rosenbluth-Putvinskii-like source 
2. the same sorce, but using the number of "fast" particles, rather than the number of runaways, to determine the source magnitude. Useful in runaway decay and/or when E<E_c, (when there are no runaways by definition).DISCLAIMER: This is not necessarily consistent, and is sensitive to the definition of a fast particle!
3. include a Chiu-Harvey-like source 
4. include a Chiu-Harvey-like source which takes thepitch-dependence of the distribution into account (derivedby Ola). Should be more correct than 3. 
5. The same as 4, but including a sink term to conserve particle energy and momentum.

.. attribute:: fastParticleDefinition

relevant when sourceMode = 2

0. Do not calculate the fast particle content, is notcompatible with sourceMode 2
1. Based on where the "tail" "begins" (where f/f_M > some threshold)Note that the shape of the source in momentum space is only recalculated if the electric field changes, and may at times notbe consistent with the beginning of the tail.
2. Based on relative speed (v/v_e > some threshold)
3. Based on absolute speed (v/c > some threshold)
4. Selects the most generous of cases 2,3 and the critical speed for runaways. Useful when you want to follow n_r closely, but still have fast particle when E is close to (or below) E_c.

.. attribute:: tailThreshold

for fastParticleDefinition = 1

.. attribute:: relativeSpeedThreshold

for fastParticleDefinition = 2 

.. attribute:: absoluteSpeedThreshold

for fastParticleDefinition = 3

.. attribute:: yCutSource

Arbitrary cutoff needed for source 5. Set to empty to useyCutSource = y_crit (currently only works with constantplasma parameters) Can also be used for source 3 and 4 (cf. sourceMode 1 vs 2)

.. attribute:: runawayRegionMode

Determines how to define the runaway region andhow to calcualate moments of f in that region.Is interpreted as 0 no mather what in sourceMode 2.

0. Isotropic (y_c=y_c(xi=1) for all xi) -- default
1. xi-dependent y_c -- more correct. Most relevant for hot-tail scenarios

.. attribute:: nPointsXiInt

resolution parameter for runawayRegionMode = 1

.. attribute:: bremsMode

include an operator for bremsstrahlung radiation reaction

0. No bremsstrahlung losses
1. Simple, continuous slowing-down force (Bakhtiari, PRL, 2005)
2. Boltzmann op., full (both low and high k, angular-dependent cross section)
3. Boltzmann op., low and high k, no angular dependence 
4. Boltzmann op., only high k, no angular dependence

