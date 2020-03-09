Reference
===================
.. class:: Reference 
    
Properties
---------------

Reference values
%%%%%%%%%%%%%%%%%%%%%%%%%

Design note: intention is that all objects containing an instance of this class also is contained in this object, think of it as pairs.
Question is if more than one instance of each class should be able to share a :class:`Reference` object.
Right now it is possible to 'hack' the construction by creating a :class:`Reference` object, and for example a :class:`TimeGrid`.
Then setting the :class:`Reference` in the :class:`TimeGrid`.
Afterwards it is possible to create a new :class:`TimeGrid` and set the already created :class:`Reference` in the new :class:`TimeGrid`.
The firstly created :class:`TimeGrid` now has a :class:`Reference` object in it with a :class:`TimeGrid` object in it which is not pointing to itself.
The 'pair' structure is then broken.
Twofixes are availible: one seperating the :class:`Reference` to a copy of itself
(new pointer) and using the copy for the old :class:`TimeGrid` or save both
:class:`TimeGrid` objects in the same :class:`Reference` object.


.. attribute:: TRef

Reference temperature

.. attribute:: nRef

Reference density

.. attribute:: nueeRef

Reference collision frequency

.. attribute:: deltaRef

Reference velocity over speed of light

.. attribute:: lnLambdaRef

Reference coulumb logarithm    


Paired Objects
%%%%%%%%%%%%%%%%%%%%

.. attribute:: physicalParams

:class:`PhysicalParams` object which also contains a pointer to this class.

.. attribute:: momentumGrid

:class:`MomentumGrid` object which also contains a pointer to this class.

.. attribute:: timeGrid

:class:`TimeGrid` object which also contains a pointer to this class.

Functions
-----------------

.. function:: Reference(TRef,nRef)

Constructor

.. function:: setPhysicalParams(this,phP)

Set the :attr:`PhysicalParams` to passed variable

.. function:: setMomentumGrid(this,mg)

Set the :attr:`MomentumGrid` to passed variable

.. function:: setTimeGrid(this,tg)

Set the :attr:`TimeGrid` to passed variable

.. function:: updateReferenceVals(this,TRef,nRef)

Update TRef, nRef and all relevant attributes in :class:`TimeGrid`, :class:`MomentumGrid` and :class:`PhysicalParams` objects in this class.
