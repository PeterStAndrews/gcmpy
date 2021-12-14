

import random
from itertools import tee, combinations
from typing import Tuple, TypeVar, Callable, List

from joint_degree import _JDS

_EDGE   = TypeVar('_EDGE', bound=Tuple[int,int])    
_NODES  = TypeVar('_NODES', bound=List[int])

##############################################
## functions to construct particular motifs ##
##############################################

def cycle_motif(nodes: _NODES) -> List[_EDGE]:
    '''Accepts a list of ints and creates a cycle of adjacent node pairs. Creates two 
    iterators and then increments one before zipping together to create tuples. Finally,
    connects the start and end of the chain toegher.
    Returns: edge list as list of tuples (int, int)'''
    nodes.sort()
    a, b = tee(nodes)
    _ = next(b, None)
    edges = list(zip(a,b))
    edges.append((nodes[0], nodes[-1]))
    return edges

def clique_motif(nodes: _NODES) -> List[_EDGE]:
    '''Accepts list of ints and creates all possible pairs which it returns as a list of tuples of int pairs'''
    nodes.sort()
    return list(combinations(nodes,2))

def diamond_motif(nodes: _NODES) -> List[_EDGE]:
    '''creates a diamond motif of 4 vertices.'''
    
    if len(nodes) < 4:
        raise("Error during motif construction - diamond_motif")

    nodes.sort()
    n0, n1, n2, n3 = nodes

    edges = cycle_motif(nodes)
    edges.append((n0,n2))
    edges.append((n1,n3))
    
    return edges

##############################################
######### network edge list object ###########
##############################################

class edge_list(object):
    '''Network represented as an edge list. The GCM class uses this 
    structure to generate networks which can then be converted to 
    other network libraries. '''
    def __init__(self):
        self._edge_list = []
        
    def add_edges_from(self, edges : List[_EDGE])->None:
        for e in edges:
            self._edge_list.append(e)

##############################################
####### datastructure to house results #######
##############################################

class output_data(object):
    '''An object to store output data from the process'''

    def __init__(self,i:int):
        self._experiment = i                # experiment index
        self._name : str = ''               # name of experiment
        self._network  : edge_list = []     # networks
    

_RESULTS   = TypeVar('_RESULTS', bound=List[output_data])  

##############################################
######## Configuration model algorithm #######
##############################################

class GCM_algorithm(object):
    """Generalised configuration model algorithm """
    _num_networks    : int                                               # number of networks to create
    _motif_sizes     : List[int]                                         # list of number of nodes in each motif
    _build_functions : List[Callable[[_NODES],List[_EDGE]]]              # list of callbacks for motif construction 
    
    def __init__(self, num_networks    : int, 
                       motif_sizes     : List[int],
                       build_functions : List[Callable[[_NODES],List[_EDGE]]] = None):
        
        self._num_networks    = num_networks 
        self._motif_sizes     = motif_sizes
        self._build_functions = build_functions

    def random_clustered_graph(self, jds : _JDS)->edge_list:
        '''Generate a random graph from a given joint degree sequence of motifs. If motif constructors 
        are not specified, ValueError is raised.
        Returns a list of edges in the graph as an edge_list object''' 

        # as default construct cliques
        if self._build_functions is None:
            raise ValueError('GCM_algorithm._build_functions is None')

        # create an empty graph
        N = len(jds)
    
        # initialise a list of lists for distinct motif topology 
        stubs = list([] for _ in range(len(jds[0])))
    
        # for each node n
        for n in range(N):
            joint_degree = jds[n]
            # for each topology ... 
            for k, k_list in enumerate(stubs):
                # append node n to the list once per unique motif of a given topology 
                for _ in range(joint_degree[k]):
                    k_list.append(n)
                
        # shuffle each stub list
        for k_list in stubs:
            random.shuffle(k_list)
        
        # create edge list object
        es = edge_list()

        # for each topology list ...
        for k, k_list in enumerate(stubs):
            motif_size = self._motif_sizes[k]
            # while there are still nodes ... 
            while k_list:
                nodes = []
                # grab required number of nodes to build the motif
                for _ in range(motif_size):
                    nodes.append(k_list.pop())

                # add the edges to the network using the builder callback
                es.add_edges_from(self._build_functions[k](nodes))    
            
        # return the graph
        return es

class ResampleJDS(GCM_algorithm):
    '''Resamples a joint degree sequence to create multiple networks.'''
    def __init__(self, num_networks    : int, 
                       motif_sizes     : List[int],
                       network_name    : str = None,
                       build_functions : List[Callable[[_NODES],List[_EDGE]]] = None):
        '''Param num_networks: integer number of networks to creates
           Param motif_sizes: list of integers for num vertices per motif
           Param network_name: string identifier/classifier for network
           Param build_functions: list of callables that take a list of node IDs and return the motif edges'''
        self._allow_rewires = True          # indicate that the networks are rewired from a single JDS sample
        self._network_name = network_name 
        super().__init__(num_networks, motif_sizes, build_functions)

    def random_clustered_graph_from_resampled_jds(self, jds : _JDS):
        '''Routine to create multiple configuration model networks from a single joint degree sequence.
        This essentially rewires a given sequence of joint degrees using the configuration model.'''
        results : _RESULTS  = []
        for i in range(self._num_networks):
            res = output_data(i)
            res._name = self._network_name
            res._network = self.random_clustered_graph(jds)
            results.append(res)
        return results
