
gcmpy: Generalised Configuration Model random graphs in Python
===================================================================

Overview
--------

``gcmpy`` is a Python library that creates random graph models according
to the generalised configuration model (GCM). Random graph models provide
an excellent framework to integrate topology with dynamics. The topology 
of a network is crucial to the outcome of a dynamical process, such as an 
epidemic, occurring over a network.

To create the networks, ``gcmpy`` creates a joint degree distribution object 
through a variety of analytical or empirical methods. Once constructed, this 
joint distribution is sampled to obtain a joint degree sequence. The joint 
sequence is then used in the GCM algorithm to create an edge list.

``gcmpy`` generates networks as edge lists and therefore can be integrated 
into any graph library such as ``networkx`` or ``iGraph``. 

Networks can be given storage tags to classify the properties for database 
look-up. 

Installation
------------

You can install ``gcmpy`` directly from PyPi using ``pip``:

.. code-block:: bash

   pip install gcmpy

The master distribution of ``gcmpy`` is hosted on GitHub. To obtain a
copy, just clone the repo:

.. code-block:: bash
    
    git clone git@github.com:PeterStAndrews/gcmpy.git
    cd gcmpy
    python setup.py install



Documentation
-------------

API documentation for ``gcmpy`` is available on `ReadTheDocs <https://peterstandrews-gcmpy.readthedocs.io/en/latest/>`_


Author and license
------------------

Copyright (c) 2021, Peter Mann <pm78@st-andrews.ac.uk>

Licensed under the `GNU General Public License v2 or later (GPLv2+) <http://www.gnu.org/licenses/gpl.html>`_.