
   
# Clique covers for gcmpy
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

from itertools import combinations
from typing import List
from random import choice

from gcmpy import edge_list
from .types import _COVER

def binom(n : int, r : int)->int:
    ''' Binomial coefficient, nCr, aka the "choose" function 
            n! / (r! * (n - r)!)
    :param n: n choose
    :param r: r
    :returns int: binomial coefficient.'''
    p = 1    
    for i in range(1, min(r, n - r) + 1):
        p *= n
        p //= i
        n -= 1
    return p

class EECC(edge_list):
    '''Creates an EECC clique cover based on the algorithm developed in (1,2) for a graph
    in an edge_list format.
    
    ----- References -----

    1) Giulio Burgio, Alex Arenas, Sergio Gómez and Joan T. Matamalas: 
       Network clique cover approximation to analyze complex contagions through group interactions, 
       Comms. Phys. (2021) in press (arXiv:2101.03618)

    2) https://github.com/giubuig/DisjointCliqueCover.jl
    '''

    def __init__(self):
        self._m0 : int = 2
        super().__init__()

    def set_max_clique_size(self, m0 : int):
        self._m0 = m0

    def limited_maximal_cliques(self)->_COVER:
        '''Calculates the maximal cliques of a graph up to order m0. Any clique of order
        m>m0 is decomposed into its sub-cliques of order m0.

        :return list of maximal cliques up to size m0'''
        indxs = []
        C = self.find_cliques()
        for c in range(len(C)):
            clique_size = len(C[c])
            if clique_size > self._m0:
                indxs.append(c)
                for nc in sorted(combinations(C[c],self._m0)): # sort by first element of sublist
                    C.append(nc)
            else:
                C[c] = sorted(C[c])

        idxs = set(range(len(C))) - set(indxs)

        C = [C[i] for i in idxs]
        C = [list(item) for item in set(tuple(row) for row in C)]

        # small difference in key for sort
        for i, c in enumerate(C):
            C[i] = sorted(c)

        return sorted(C, key=lambda x: (-len(x), x[0], x[1]) if len(x) > 1  else (-len(x), x[0], 0))


    def compute_scores(self, C: _COVER,
                             EC : _COVER,
                             ord : List[int],
                             r : List[float],
                             indexes : List[int])->None:
        '''Scores the cliques in list C in the network G and includes those of score zero in the 
        cover EC. The higher the score
        
        :param C: set of maximal clique
        :param EC: EECC cover
        :param ord: List of cliques' order
        :param r: List of cliques' score
        :param indexes: List of indexes'''

        num_cliques = len(C)
        for c in range(num_cliques):
            
            order = len(C[c])   # num vertives in clique
            C[c] = sorted(C[c]) # sort vertices
            ord[c] = order      # set size
            
            if order > 2:

                size = binom(order,2) # number of edges in clique of size order

                for i in range(order): # for each vertex pair in clique (i.e. edges)
                    for j in range(i+1,order):

                        f = 0
                        n = 0
                        # iterate all cliques and see if *this* edge (i,j)
                        # is a part of other cliques or not
                        while n <= num_cliques-1:
                            if n != c and set([C[c][i],C[c][j]]).issubset(C[n]):
                                f = 1
                                break
                            else:
                                n+=1

                        if f == 1:
                            # increase the score for the overlapping edge
                            r[c] += 1.0/size

            # if score is zero, add the clique to the cover      
            if r[c] == 0:
                EC.append(C[c])
                indexes.append(c)

    def get_EECC(self)->_COVER:
        '''Calculate the edge-disjoint edge clique cover (EECC) of a graph G considering 
        cliques of order up to m0, according to the heuristic proposed in reference (1).
        
        :param G: graph
        :param m0: int for maximum order of the cliques to consider

        :return cover: A list containing the EECC'''

        C = self.limited_maximal_cliques()
        num_cliques : int = len(C)
        ord : List[int] = [0] * num_cliques
        r : List[float] = [0.0] * num_cliques
        EC : _COVER = []

        indexes_score0 = []
        self.compute_scores(C,EC,ord,r,indexes_score0)

        len_ = len(C)
        idxs = set(range(len_)) - set(indexes_score0)

        Ctemp = []
        ordtemp = []
        rtemp = []
        for i in idxs:
            Ctemp.append(C[i])
            ordtemp.append(ord[i])
            rtemp.append(r[i])

        C = Ctemp
        ord = ordtemp
        r = rtemp
        for c in range(len(EC)):
            order = len(EC[c])
            for i in range(order):
                for j in range(i+1,order):
                    self.remove_edge(EC[c][i],EC[c][j])

        while self.has_edges():
            # find (one of) the largest clique(s) with the minimum score, include it in EECC

            # get smallest score
            min_r : float = min(r)

            # get indices of all smallest scoring cliques 
            min_r_set = [idx for idx, element in enumerate(r) if element == min_r]
            
            # find size of largest clique of the smallest scoring cliques
            max_ord : int = -1
            for idx in min_r_set:
                if ord[idx] > max_ord:
                    max_ord = ord[idx]

            # construct a list of lowest-scoring, largest sized cliques
            indexes_to_sample = [idx for idx in min_r_set if ord[idx] == max_ord]

            # stochastic part of algorithm 
            idx = choice(indexes_to_sample)

            # the selected clique to include in the EECC
            cli = C[idx]
            EC.append(cli)

            # remove the clique from the graph
            for i in range(max_ord):
                for j in range(i+1,max_ord):
                    # assumes edges are ordered i < j
                    self.remove_edge(cli[i],cli[j])

            C = self.limited_maximal_cliques()

            # remove isolated vertices
            C = [c for c in C if len(c) > 1]
            ord = [0] * len(C)
            r = [0] * len(C)
            indexes_score0 = []

            self.compute_scores(C, EC, ord, r, indexes_score0)
            
            # Removing those cliques with score zero from C, r and G_
            len_ = len(C)
            idxs = set(range(len_)) - set(indexes_score0)

            Ctemp = []
            ordtemp = []
            rtemp = []
            for i in idxs:
                Ctemp.append(C[i])
                ordtemp.append(ord[i])
                rtemp.append(r[i])


            C = Ctemp
            ord = ordtemp
            r = rtemp
            for c in range(len(EC)):
                order = len(EC[c])

                for i in range(order):
                    for j in range(i+1,order):
                        self.remove_edge(EC[c][i],EC[c][j])

        return sorted(EC, key=lambda x: (-len(x), x[0], x[1]))


class MPCC(edge_list):
    '''Creates an MPCC clique cover. '''

    def __init__(self):
        pass    
    
