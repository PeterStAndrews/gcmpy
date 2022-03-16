

import random
import networkx as nx
from iteration_utilities import grouper
from itertools import chain,repeat,starmap 

from gcmpy.network.network import Network

class GCMAlgorithm:
    """
    Generalised configuration model algorithm.
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

    def random_clustered_graph(self, jds: list) -> Network:
        """
        Generate a random graph from a given joint degree sequence of motifs. If motif constructors 
        are not specified, ValueError is raised.
        :param jds: joint degree sequence 
        :returns: a list a networkx graph object
        """
        stubs = [list(chain.from_iterable(starmap(repeat,r)))
           for r in map(enumerate,zip(*jds))  ]
    
        # shuffle each stub list
        for k_list in stubs:
            random.shuffle(k_list)

        # create graph object
        model = Network()

        model._G.add_nodes_from([i for i in range(len(jds))])

        # give vertices joint degree attribute
        joint_degree_dict = {n : [0]*len(self._motif_sizes) for n in model._G.nodes()}
        nx.set_node_attributes(model._G, joint_degree_dict, 'joint_degree')

        #for each topology list ...
        for k, k_list in enumerate(stubs):
            # iterate the degree list
            for nodes in grouper(k_list,self._motif_sizes[k]):
                # add the edges to the network using the builder callback
                es = self._build_functions[k](list(nodes))
                model._G.add_edges_from(es)

                # update the joint degrees of each vertex
                for n in nodes:
                    model._G.nodes[n]['joint_degree'][k] += 1

                # give topology name to edges
                if self._edge_names:
                    for e in es:
                        model._G.edges[e]['topology'] = self._edge_names[k]     
            
        # return the graph model
        return model