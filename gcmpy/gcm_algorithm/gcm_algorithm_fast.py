

import random
from iteration_utilities import grouper
from itertools import chain,repeat,starmap

class FastGCMAlgorithm:
    """
    Generalised configuration model algorithm. Fast version doesn't use 
    networkx or store vertex and edge data.
    :param num_networks: the number of networks to create
    :param motif_sizes: list of ints that indicate the number of nodes in each motif
    :param build_functions: callbacks that accept list of nodes and return edges
    """
    _num_networks: int                      # number of networks to create
    _motif_sizes: list                      # list of number of nodes in each motif
    _build_functions: list                  # list of callbacks for motif construction 
    _edge_names: list                       # list of edge topology names
    
    def __init__(self, motif_sizes: list,
                       build_functions: list,
                       edge_names: list = None):
        
        self._motif_sizes = motif_sizes
        self._build_functions = build_functions
        self._edge_names = edge_names

    def random_clustered_graph(self, jds: list) -> list:
        """
        Generate a random graph from a given joint degree sequence of motifs. If motif constructors 
        are not specified, ValueError is raised.
        :param jds list: joint degree sequence 
        :returns list: a list of edges
        """
        stubs = [list(chain.from_iterable(starmap(repeat,r)))
           for r in map(enumerate,zip(*jds))  ]
    
        # shuffle each stub list
        for k_list in stubs:
            random.shuffle(k_list)

        # create list for edges
        edge_list = []

        #for each topology list ...
        for k, k_list in enumerate(stubs):
            # iterate the degree list
            for nodes in grouper(k_list,self._motif_sizes[k]):
                # add the edges to the network using the builder callback
                es = self._build_functions[k](list(nodes))
                edge_list.extend(es)

        # return the graph model
        return edge_list