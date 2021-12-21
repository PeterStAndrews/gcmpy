.. gcmpy documentation master file, created by
   sphinx-quickstart on Sat Dec 18 23:31:24 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to gcmpy: a toolkit for random graph creation for network science
===========================================================================

``gcmpy`` aims to provide a common framework for the scalable and
efficient creation of random graphs in python that are constructed 
according to the generalised configuration model algorithm.

What are configuration model networks?
---------------------------------------

Networks are collections of vertices connected together by edges. A mathematician 
might call a network a *graph*, and we use the two terms indiscriminately. Networks provide 
an excellent opportunity to bring both dynamics and topology into a single model. In this context
questions on how topology influences the dynamics are readily answered by a network model. For instance,
in a *contact network*, the vertices of the graph represent people, and the edges represent those 
contacts that are sufficient to transmit disease. One could then inject a dynamical process over 
this network, such as an epidemic disease that transmits from contact to contact. Varying the structure
of the contact network, perhaps through mitigation strategies such as social distancing, might change
how the epidemic behaves at the population level. 

Networks can be collected from empirical data, or they can be created according to a given algorithm.
Empirical networks are only a single graph; whilst, synthetic algorithms can be used repeatedly to 
create different graphs. In general, graph construction algorithms are stochastic by nature, meaning 
that they have randomisation in them. This leads to a different final network for a given input dataset. 
As such, a *family* of random graphs can be created with equivalent statistical properties. Upon each 
application of the algorithm, a single member of the random graph family is generated.

One such random graph algorithm is the *configuration model* (CM). In the CM, each vertex in a set 
of unconnected vertices is assigned an integer, :math:`k`, that represent it's number of edges - called its degree.
The vertices are then entered into a list once for each of its edges, such that, a vertex with :math:`k` edges 
is present exactly :math:`k` times in the list. Once all vertices have been entered, the list is shuffled
and vertex pairs are drawn without replacement from the list. The vertices are then connected together. Once 
all pairs have been drawn from the list, the construction algorithm is completed. Theoretically, the same vertex
can be drawn from the list twice, and therefore, self-loops can be created; however, for large networks, this
is rare.

What are generalised configuration model networks?
---------------------------------------------------

In the CM, pairs of vertices are connected together at random. This leads to networks that are 
*locally tree-like* and do not contain short range cycles. In the contact network example above
this indicates that an individuals friends do not have contact with one another, which is highly 
unrealistic. These types of relationships often have a very large influence on the properties of 
the network!

The *generalised configuration model* (GCM) allows higher-order connections, beyond simple edges, to be added to the graph whilst
still generating random graphs. An edge is a pair-wise interaction. The next interaction we can consider is 
an interaction between 3 vertices, such as a triangle. Let us consider a GCM network that contains both ordinary edges 
and triangles. To begin, each vertex in an unconnected set of vertices is assigned two integers: one representing
its number of ordinary edges, :math:`s`, and the second representing the number of distinct triangles, :math:`t`, it belongs to. Two lists are 
now created and each vertex is entered once for each of its ordinary edges to the first list, and once for each of its 
triangles to the second. Once all vertices have been added to the lists they are sboth huffled. In a similar 
manner to the CM algorithm, pairs of vertices are drawn from the ordinary edge list and connected together until it is empty; whilst triples
of vertices are drawn from the triangles list and connected together in a triangle until it too is empty. Upon completion
we now have a random network that contains ordinary edges and triangles. 

The choice of motif can certainly be extended beyond ordinary edges and triangles. Indeed, any graph motif can be used.
This allows a great freedom in the kinds of networks that we can create and the types of natural phenomena we can model!

What is ``gcmpy``?
---------------------

``gcmpy`` is a pure Python simulation framework that facilitates the synthesis of GCM networks from a variety of 
input formats. It aims to provide the common simulation approaches used in
the scientific literature, together with a small set of "common
networks" that can form the basis for experimentation.

Features
--------

* Compatible with Python 3.8 and later

* Fully compatible with ``jupyter`` notebooks and labs

* Annotated with ``typing`` type annotations

* Supports a variety of analytical and empircal entry points for network construction

* Allows empirical networks to be covered and synthetic data to be generated

* Supports serialisation for database storage and tags for retrieval

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Contents:

   tutorial
