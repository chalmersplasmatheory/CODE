PhysicalParams:
====================

Class containing the physical parameters for the run.
The idea is that the parameters you get out will always be interpolated to a :class:`TimeGrid`.
This is why you only are allowed to change parameters through function calls as this allows for consistent updates with :class:`TimeGrid`.
The reason for this is that if something changes, it is best if only you need to look at one place.
So let us say that this would be done in the relevant implementation of :class:`Solver` class instead, then if something needs special treatment in
this class, you would need to update all implementations of :class:`Solver`.
However, as it is implemented now it is hard to change a variables name (as this is fetched by this class property) but you always now that it is 
consistent with the :class:`TimeGrid`.
Therefore, there are two sets for each supplied physical parameter, the raw supplied parameter from which the interpolation is done and the 
interpolated parameter. Example :attr:`rawT` which is the raw usersuppplied temperature and :attr:`T` which is its interpolated value.

.. class:: PhysicalParams

Properties
----------------

Physical quantities:
%%%%%%%%%%%%%%%%%%%%%

.. attribute:: T 

Plasma temperature (eV), given in timesteps

.. attribute:: n 

Plasma density (m^{-3}), given in timesteps

.. attribute:: Z 

Plasma effective charge, given in timesteps

.. attribute:: E 

Electric field (V/m), given in timesteps

.. attribute:: B 

(optional) Magnetic field (T), only used for calculating
synchrotron radiation reaction. If a magnetic field strength is
provided, a term describing momentum loss due to synchrotron
emission is included in the kinetic equation. If B is empty it is
not included.

.. attribute:: species 

a struct used with the species properties so that species has the fields nj, Z0NetCharge, ZAtomicNumber, times:

* nj(jspecies, times) -  matrix containing number density of each species as a function of time (m :sup:`-3`)
* ZAtomicNumber -  atomic number of each species. row vector, not a function of time
* Z0NetCharge - net charge for each species. row vector, not a function of time
* times   time steps, same as the timessteps vector in :attr:`timeGrid`.
* NB species cannot be set from SetPhysicalParameters. 
  NB overrides specified :attr:`Z` and :attr:`n` since nj and Z0NetCharge specifies the effective charge. 
  (Z is the fully screened value) NB only used if useScreening ~= 0

.. attribute:: rawspecies

Usersupplied species struct (non time interpolated)

a struct used with the species properties so that species has the fields nj, Z0NetCharge, ZAtomicNumber, times:

* nj(jspecies, times) -  matrix containing number density of each species as a function of time (m :sup:`-3`)
* ZAtomicNumber -  atomic number of each species. row vector, not a function of time
* Z0NetCharge - net charge for each species. row vector, not a function of time
* times   time steps, where the ionization and densities are set
* NB species cannot be set from SetPhysicalParameters. 
  NB overrides specified :attr:`Z` and :attr:`n` since nj and Z0NetCharge specifies the effective charge. 
  (Z is the fully screened value) NB only used if useScreening ~= 0

.. attribute:: rawtT 

times where rawT is defined, user supplied but to normalized time units (normalized)

.. attribute:: rawtn 

times where rawn is defined, user supplied but to normalized time units (normalized)

.. attribute:: rawtZ 

times where rawZ is defined, user supplied but to normalized time units (normalized)

.. attribute:: rawtE 

times where rawE is defined, user supplied but to normalized time units (normalized)

.. attribute:: neTotalOverneFree 

vector containing fraction density of total
electrons compared to free electrons at all timesteps.
Used when specicies is used (for impurities and boltzmann
operator, and bound electrons are present). If no species are
used or no source is used, this is set to 1.

.. attribute:: nuees 

Collisional frequency of electron electron collisions

.. attribute:: deltas 

thermal speed over speed of light

.. attribute:: EOverEc 

Electric field over the Critical Electrical field

.. attribute:: EOverED 

Electric field over the Dreicer Electrical field

.. attribute:: lnLambdas 

.. math::
        \ln \Lambda = 14.9 - 
                0.5 \ln (\frac{n}{10^{20} \text{ m}^{-3}}) +
                \ln (\frac{T}{10^3 \text{ eV}})


Coloumb logarithm

Miscellaneous:
%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: interpolationMethod 

Interpolation method to interpolate between
given values and values at timestepping times.
Used in interp1 (default previous) NOTE: SINCE
NEAREST IS DEFAULT, T WILL DROP AT HALF THE
TIME YOU SPECIFY WHERE T SHOULD DROP

Time Grid:
%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: timeGrid 

:class:`TimeGrid` object. 
It is at these timesteps the parameters are set.

Reference values:
%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: reference

:class:`Reference` object.

.. attribute:: VERSION 

Version of this object

Functions:
--------------------

.. function:: this = PhysicalParams(reference,timeGrid,varargin)

Set Physical properties:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. function:: setParams(this,varargin)

Set params in standard matlab syntax, *e.i.* 'name' followed by value.  
Does not support update of reference or TimeGrid.

.. function:: setPhysicalParameters(o,T,n,Z,E,varargin)

Usage:

* setPhysicalParameters(T,n,Z,E)
* setPhysicalParameters(T,n,Z,E,B)
* setPhysicalParameters(T,n,Z,E,B,tT,tn,tZ,tE)

t* should be in normalized units.
If t* is in seconds, simply multiply by :attr:`reference`.nueeRef.

.. function:: setspecies(this,species)

Sets the species struct

.. function:: setT(this,T,tT)

Sets the temperature. 
tT is the time at which the temperature is defined.
tT is optional if T is scalar.
tT should be in normalized units (simply if tn is in seconds, mulitply by :attr:`reference`.nueeRef)

.. function:: setn(this,n,tn)

Sets the density (m^-3). 
tn is the time at which the density is defined.
tn is optional if n is scalar.
tn should be in normalized units (simply if tn is in seconds, mulitply by :attr:`reference`.nueeRef)

.. function:: setZ(this,Z,tZ)

Sets the charge of plasma. 
tZ is the time at which the charge is defined.
tZ is optional if Z is scalar.
tZ should be in normalized units (simply if tn is in seconds, mulitply by :attr:`reference`.nueeRef)

.. function:: setE(this,E,tE)

Sets the Electric field (V/m)
tE is the time at which the electric field is defined.
tE is optional if E is scalar.
tE should be in normalized units (simply if tn is in seconds, mulitply by :attr:`reference`.nueeRef)

.. function:: setB(this,B)

Sets the magnetic field (Teslas).

Calculate Dependent Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. function:: calcDepParams(this)

Misc.
%%%%%%%%%%

.. function:: interpolatePhysicalParams(this)

.. function:: updateRawTimes(this, nueeRefNew)

PHYSICALPARAMS DO NOT USE THIS FUNCTION, IT IS ONLY FOR
REFERENCE CLASS. NO SOLUTION FOUND WHERE ONLY REFERENCE CLASS
CAN GET THE RAW TIMES FOUND WHERE YOU ALSO CAN USE THE
FUNCTION. DONT USE THIS FUNCTION UNLESS IN REFERENCE CLASS
