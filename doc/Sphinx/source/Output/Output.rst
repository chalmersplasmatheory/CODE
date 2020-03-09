Output
============
.. class:: Output

Properties
----------------

Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: saveDist 

switch if distribution is saved, default true

Saved Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. attribute:: distributions

Cell with :class:`Distributions` at times. 

.. attribute:: yc

Cell with yc at times.

.. attribute:: nr

Cell with nr at times.

.. attribute:: averageEnergyMeV

Cell with averageEnergyMeV at times.

.. attribute:: averageFastEnergyMeV

Cell with averageFastEnergyMeV at times.

.. attribute:: averageREEnergyMeV

Cell with averageREEnergyMeV at times.

.. attribute:: currentDensity

Cell with currentDensity at times.

.. attribute:: currentDensityRE

Cell with currentDensityRE at times.

.. attribute:: fracRE

Cell with fracRE at times.

.. attribute:: fracRECurrent

Cell with fracRECurrent at times.

.. attribute:: fracREEnergy

Cell with fracREEnergy at times.

.. attribute:: growthRate

Cell with growthRate at times.

.. attribute:: growthRatePerSecond

Cell with growthRatePerSecond at times.

.. attribute:: times

Cell with times where other values are saved.

.. attribute:: totalEnergyMeV

Cell with totalEnergyMeV at times.

.. attribute:: totalEnergyJ

Cell with totalEnergyJ at times.

.. attribute:: dJdtSI

Cell with dJdtSI at times.

.. attribute:: dJ_rundtSI

Cell with dJ_rundtSI at times.

.. attribute:: T

Cell with T at times.

.. attribute:: n

Cell with n at times.

.. attribute:: Z

Cell with Z at times.

.. attribute:: E

Cell with E at times.

.. attribute:: B

Cell with B at times.

.. attribute:: species

Cell with species at times.

.. attribute:: neTotalOverneFree

Cell with neTotalOverneFree at times.

.. attribute:: nuees

Cell with nuees at times.

.. attribute:: deltas

Cell with deltas at times.

.. attribute:: EOverEc

Cell with EOverEc at times.

.. attribute:: EOverED

Cell with EOverED at times.

.. attribute:: EHats

Cell with EHats at times.

.. attribute:: BHatRef

Cell with BHatRef at times.

.. attribute:: nueeBars

Cell with nueeBars at times.

.. attribute:: nBars

Cell with nBars at times.

.. attribute:: veBars

Cell with veBars at times.

.. attribute:: veBars

Cell with veBars at times.

.. attribute:: veBars

Cell with veBars at times.

.. attribute:: lnLambdas

Cell with lnLambdas at times.


Functions
-------------

.. function:: this = Output()

.. function:: save(this, timeIndex, state, f, fbefore)

Add another entry in all saved properties.

.. function:: addDistribution(this,f,state,saveIndex)

Adds distrubtion function to saved cells, where values of distribution corresponds to saveIndex of state.physicalParams.[PARAM](saveIndex)

.. function:: [f, times] = getDistributions(this)

Returns distributions and their corresponding times into matlab vector format.
Each coloumn is a distribution function at a specific time.

.. function:: [p, times] = getMomentumVectors(this)

Returns the momentum vectors used for the distributions functions returned in :func:`getDistributions` and their respective times
