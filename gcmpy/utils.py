# Utilities for gcmpy
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

import networkx as nx
from .types import _EDGES, _EDGE

class network(object):
    '''Network represented as an edge list. The GCM class uses this 
    structure to generate networks which can then be converted to 
    other network libraries.'''

    def __init__(self):
        self._G = nx.Graph()
        
    def add_edge(self, e : _EDGE):
        '''Adds an edge to the networkx'''
        self._G.add_edge(*e)

    def add_edges_from(self, edges : _EDGES)->None:
        '''Adds edges from list of tuples (int,int) to the edge list.
        :param edges: list of tuples of ints.'''
        self._G.add_edges_from(edges)

    def find_cliques(self):
        '''Returns all maximal cliques'''
        return list(nx.find_cliques(self._G))

    def remove_edge(self, i : int, j : int)->None:
        '''Removes edge (i,j) from G. With nx, will also remove (j,i).'''
        try:
            self._G.remove_edge(i,j)
        except nx.NetworkXError:
            return
        
    def has_edges(self)->bool:
        '''True if graph has edges remaining'''
        return len(self._G.edges()) > 0