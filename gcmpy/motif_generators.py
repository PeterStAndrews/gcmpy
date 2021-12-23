# Motif generators for gcmpy
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

from itertools import tee, combinations
from typing import List
from .types import _NODES, _EDGE

def cycle_motif(nodes: _NODES) -> List[_EDGE]:
    '''Accepts a list of ints and creates a cycle of adjacent node pairs. Creates two 
    iterators and then increments one before zipping together to create tuples. Finally,
    connects the start and end of the chain toegher.
    
    :param nodes: list of nodes (int)
    :returns: edge list as list of tuples (int, int)'''
    #nodes.sort()
    a, b = tee(nodes)
    _ = next(b, None)
    edges = list(zip(a,b))
    edges.append((nodes[0], nodes[-1]))
    return edges

def clique_motif(nodes: _NODES) -> List[_EDGE]:
    '''Accepts list of ints and creates all possible pairs which it returns as a list of tuples of int pairs
    
    :param nodes: list of nodes (int)
    :returns: edge list as list of tuples (int, int)'''
    #nodes.sort()
    return list(combinations(nodes,2))

def diamond_motif(nodes: _NODES) -> List[_EDGE]:
    '''Accepts list of ints and creates edge pairs for a diamond motif of 4 vertices.
    
    :param nodes: list of nodes (int)
    :returns: edge list as list of tuples (int, int)'''
    if len(nodes) < 4:
        raise("Error during motif construction - diamond_motif")

    #nodes.sort()
    n0, n1, n2, n3 = nodes

    edges = cycle_motif(nodes)
    edges.append((n0,n2))
    edges.append((n1,n3))
    
    return edges
