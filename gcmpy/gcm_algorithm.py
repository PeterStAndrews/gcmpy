# GCM algorithm for gcmpy
#
# Copyright (C) 2021 Peter Mann
#
# This file is part of gcmpy, generalised configuration model networks in Python.
#
# gcmpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# gcmpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gcmpy. If not, see <http://www.gnu.org/licenses/gpl.html>.

import random
from typing import List, Callable
from .types import _EDGES, _NODES, _JDS
from .utils import edge_list, output_data, results

class GCM_algorithm(object):
    """Generalised configuration model algorithm.
    
    :param num_networks: the number of networks to create
    :param motif_sizes: list of ints that indicate the number of nodes in each motif
    :param build_functions: callbacks that accept list of nodes and return edges"""
    _num_networks    : int                                               # number of networks to create
    _motif_sizes     : List[int]                                         # list of number of nodes in each motif
    _build_functions : List[Callable[[_NODES],_EDGES]]              # list of callbacks for motif construction 
    
    def __init__(self, motif_sizes     : List[int],
                       build_functions : List[Callable[[_NODES],_EDGES]]):
        
        self._motif_sizes     = motif_sizes
        self._build_functions = build_functions

    def random_clustered_graph(self, jds : _JDS)->edge_list:
        '''Generate a random graph from a given joint degree sequence of motifs. If motif constructors 
        are not specified, ValueError is raised.
        
        :param jds: joint degree sequence 
        :returns: a list of edges in the graph as an edge_list object''' 

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
    '''Resamples a joint degree sequence to create multiple networks.
    
    :param num_networks: the number of networks to create
    :param motif_sizes: list of ints that indicate the number of nodes in each motif
    :param network_name: string identifier/classifier for network
    :param build_functions: callbacks that accept list of nodes and return edges'''
    def __init__(self, num_networks    : int, 
                       motif_sizes     : List[int],
                       build_functions : List[Callable[[_NODES],_EDGES]],
                       network_name    : str = None,):
        self._num_networks = num_networks
        self._allow_rewires = True          # indicate that the networks are rewired from a single JDS sample
        self._network_name = network_name 
        super().__init__(motif_sizes, build_functions)

    def random_clustered_graph_from_resampled_jds(self, jds : _JDS)->results:
        '''Routine to create multiple configuration model networks from a single joint degree sequence.
        This essentially rewires a given sequence of joint degrees using the configuration model.
        
        :param jds: joint degree sequence 
        :returns results: data of constructed networks'''
        res : results = []
        for i in range(self._num_networks):
            res_ = output_data(i)
            res_._name = self._network_name
            res_._network = self.random_clustered_graph(jds)
            res.append(res_)
        return res
