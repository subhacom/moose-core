.. _intermediate-python:

********************************************
More on python scripting for MOOSE
********************************************

.. _load-model:

Loading predefined models
=========================

.. _searching-with-wildcardfind:

Searching the element tree
==========================

In the very first example in the beginners' tutorial
(:ref:`importing-moose`) you saw ``moose.wildcardFind`` function::

  >>> for x in moose.wildcardFind( '/model/#graphs/conc#/#' ):
  ...     pylab.plot( x.vector, label=x.name )

Whether you loaded an existing model or set up your own in Python,
typing out the path of every object in the model or keeping track of
every object that you create in Python variables become tedious with
larger models. ``moose.wildcardFind`` is a powerful function that
allows you to search the moose element tree for objects that meet some
criterion. MOOSE uses `'#'` as the wildcard character in paths. This
is similar to `'*'` used for file searches in most operating systems.

`wildcardFind` takes a string argument which represents a path with a wildcard::

  >>> moose.wildcardFind('/#')
  
This will return a list of elements under `root`, i.e., its children. A single `'#'`
character after an element's path matches all child elements.

A `'##'` triggers recursive search under the path::

  >>> moose.wildcardFind('/##')

This will show a long list of elements including the elements we
created under `/model`.

You can match parts of the name of an element in the path::

  >>> moose.wildcardFind('/##/stim#')
  [<moose.PulseGen id=489 dataIndex=0 path=/model[0]/stimulus[0]>]


You can restrict the match by adding conditions after the wildcard::

  >>> moose.wildcardFind('/##[TYPE=Compartment]')
  [<moose.Compartment id=487 dataIndex=0 path=/model[0]/soma[0]>]

Here `[TYPE=Compartment]` tells `wildcardFind` to select only those
elements which are of ``moose.Compartment`` class.

You can also set the condition to compare the field value to a specific value::

  >>> moose.wildcardFind('/##[FIELD(Vm)>-0.07]')
  [<moose.Compartment id=541 dataIndex=0 path=/model[0]/soma[0]>]


In this case we have just one compartment. In a multicompartmental
model, or networks of multicompartmental neurons, this kind of search
is handy.


.. _moose-messagaes:

Messaging in MOOSE: How elements pass information
=================================================

moose.showmsg()


.. _neighbor-elements:

Finding neighbors
=================
moose.neighbors()
element.neighbors[field]
