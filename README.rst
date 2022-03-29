gcmpy: Generalised Configuration Model random graphs in Python
===================================================================

.. image:: https://badge.fury.io/py/gcmpy.svg
    :target: https://badge.fury.io/py/gcmpy

.. image:: https://readthedocs.org/projects/peterstandrews-gcmpy/badge/?version=latest
    :target: https://peterstandrews-gcmpy.readthedocs.io/en/latest/?badge=latest
      
.. image:: https://github.com/PeterStAndrews/gcmpy/actions/workflows/ci.yml/badge.svg
     :target: hhttps://github.com/PeterStAndrews/gcmpy/actions/workflows/ci.yml

.. image:: https://static.pepy.tech/personalized-badge/gcmpy?period=total&units=international_system&left_color=black&right_color=green&left_text=Downloads
    :target: https://pepy.tech/project/gcmpy

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
sequence is then used in the GCM algorithm to create a networkx graph. It can 
also be used to create an edge list directly, which is significantly faster.

There is also a tools library for obtaining useful quantities from the network
as well as converting a joint degree distribution into excess joint degree 
distributions and vice versa, for example. We also provide an MCMC rewiring algorithm 
to stochastically rewire a synthetic network's correlations structure to a  
target joint excess joint degree distributions. 

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


The unit tests can be discovered from the project root using 

.. code-block:: bash

    python3 -m unittest discover -v -s test/ -p 'test_*.py'

Documentation
-------------

API documentation for ``gcmpy`` is available on `ReadTheDocs <https://peterstandrews-gcmpy.readthedocs.io/en/latest/>`_


Author and license
------------------

Copyright (c) 2021, Peter Mann <pm78@st-andrews.ac.uk>

Licensed under the `GNU General Public License v2 or later (GPLv2+) <http://www.gnu.org/licenses/gpl.html>`_.