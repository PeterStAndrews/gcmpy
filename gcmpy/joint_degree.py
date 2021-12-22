# construct joint degree distributions for gcmpy
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

import ast
import random
from itertools import product
from collections import Counter

import numpy as np

from typing import Tuple, Callable, List
from .types import _JDD, _JDS, _JOINT_DEGREE, _COVER

class JDD_Interface(object):
    '''Joint degree distribution interface. Subclasses will define how self._jdd
    is created, and all implementations should do so upon construction.
    
    :param motif_sizes: list of integers for number of nodes in each motif'''

    def __init__( self, motif_sizes ):
        self._jdd : _JDD = {}
        self._motif_sizes : List[int] = motif_sizes

    def _sample_JDS( self, N : int )->_JDS:
        '''Samples N joint motif sequence from self._jdd object. Intended to 
        be a private method, only to be called by sample_JDS.
        
        :param N: number of samples
        :return jds: joint degree sequence'''
        
        keys = np.array(list(self._jdd.keys()))
        values = np.array(list(self._jdd.values()))
        choice_list_indices = np.random.choice(len(keys), N, replace=True, p=values)
        return keys[choice_list_indices] 

    def sample_JDS( self, N : int )->_JDS:  
        '''Method to ensure the sum of the motifs in each dimension mod motif_sizes[i] (which is the 
        number of vertices in each motif) is zero by adding motifs to the motif list. 
        This ensures the JDS is graphic.
        
        :param N: number of samples
        :returns jds: joint degree sequnce'''
        
        # raw joint degree sequence
        jds = list(map(ast.literal_eval,self._sample_JDS(N)))
        
        # ensure that the sum of the numbers is divisible by the motif size
        ntops = list(map(sum, zip(*jds)))
        for i, ntop in enumerate(ntops):
            if ntop % self._motif_sizes[i] != 0:
                # if not, round up to add one more motif to the network
                for j in range(self._motif_sizes[i] - ntop % self._motif_sizes[i]):
                    j = random.randrange(0,len(jds))
                    t = list(jds[j])
                    t[i] += 1
                    jds[j]=t
        return jds

    def normalise_jdd(self):
        '''Normalises self._jdd probability distribution.'''
        summation = sum(self._jdd.values())
        for key in self._jdd:
            self._jdd[key] /= summation

    def convert_jds_to_jdd(self, jds : _JDS):
        '''Convert a joint degree sequence to self._jdd.
        
        :param jds: joint degree sequence'''
        n_samples = len(jds)
        self._jdd = dict((k,v/n_samples) for k, v in Counter(str(tuple(e)) for e in list(jds)).items())
 
class JDD_manual(JDD_Interface):
    '''Specify self._jdd by hand.
    :param jdd: joint degree distribution''' 
    
    def __init__(self, jdd : _JDD, motif_sizes : List[int])->None:
        super().__init__(motif_sizes)
        self._jdd = jdd
        

class JDD_empirical_data(JDD_Interface):
    '''An empirical joint degree sequence is used to create self._jdd.
    :param jds: joint degree sequence'''

    def __init__(self, jds : _JDS, motif_sizes : List[int])->None:
        super().__init__(motif_sizes)
        self.convert_jds_to_jdd(jds)
        
class JDD_joint_function(JDD_Interface):
    '''Multivariate function to evaluate the probability of given joint degree from an analytical source.
    Note, callable self._fp must accept a joint degree tuple and return a float.
    
    :param fp: callback
    :param motif_sizes: list of ints for number of vertices in each motif
    :param hi_lo_degree_bounds: list of tuples (int,int) for kmin,kmax per topology
    :param use_sampling: bool to use sampling or direct approach
    :param n_samples: number of samples if not direct'''
    def __init__(self, fp : Callable,
                       motif_sizes : List[int],
                       hi_lo_degree_bounds : Tuple[int,int],
                       use_sampling : bool = False,
                       n_samples : int = 1e5
                       ):
        super().__init__(motif_sizes)
        self._fp = fp
        self._hi_lo_degree_bounds = hi_lo_degree_bounds  
        
        if not use_sampling:
            self._create_jdd_directly()
        else:
            self._n_samples = n_samples
            self._create_jdd_by_sampling()

    def _create_jdd_by_sampling(self)->None:
        '''Sampling algorithm for joint degree distribution function not defined. Should 
        be overriden for each implementation.'''
        raise(NotImplementedError)

    def _create_jdd_directly(self)->None:
        '''Evaluates probability directly by generating all possible joint degrees.'''
        # build list of lists of possible degrees in each dimension
        ks = [list(range(kmin,kmax+1)) for kmin,kmax in self._hi_lo_degree_bounds]

        # iterate all joint degrees and evaluate the joint degree
        for jd in list(product(*ks)):
            self._jdd[str(jd)] = self._fp(jd)

class JDD_marginals(JDD_Interface):
    '''Merge uncorrelated marginals in each topology from analytical data together to create self._jdd. 
    If using a direct method, all possible joint degree tuples are evaluated; however, for large varience 
    in the allowed degrees this method is slow. Instead, we can choose to sample the analytical functions
    by setting use_sampling which draws n_samples weighted samples from each marginal function.
    
    :param arr_fp: array of callbacks
    :param motif_sizes: list of ints for number of vertices in each motif
    :param hi_lo_degree_bounds: list of tuples (int,int) for kmin,kmax per topology
    :param use_sampling: bool to use sampling or direct approach
    :param n_samples: number of samples if not direct'''
    def __init__(self, arr_fp : List[Callable],
                       motif_sizes : List[int],
                       hi_lo_degree_bounds : List[Tuple[int,int]],
                       use_sampling : bool = False,
                       n_samples : int = 1e5
                       ):
        super().__init__(motif_sizes)
        self._arr_fp = arr_fp
        self._hi_lo_degree_bounds = hi_lo_degree_bounds

        if not use_sampling:
            self._create_jdd_directly()
        else:
            self._n_samples = n_samples
            self._create_jdd_by_sampling()

    def _generate_all_joint_degrees(self)->List[_JOINT_DEGREE] :
        '''Generate all possible joint degrees from a range of min/max degree of each motif.
        :return jd: list of joint degrees'''
        ks = []
        for kmin, kmax in self._hi_lo_degree_bounds:
            ks.append([k for k in range(kmin,kmax)])
        return list(product(*ks))
    
    def _evaluate_prob_of_joint_degree(self, joint_degree : List[_JOINT_DEGREE])->float:
        '''Evaluate the joint probability of a joint degree using function callbacks.
        
        :param joint_degree: list of joint degreees
        :returns prod: probability of joint degree'''
        prod : float = 1.0
        for i,deg in enumerate(joint_degree):
            prod *= self._arr_fp[i](deg)
        return prod

    def _create_jdd_directly(self)->None:
        '''create self._jdd by enumerating probability of all possible joint degrees from 
        analytical functions.'''

        self._jdd = dict((str(key), 0.0) for key in self._generate_all_joint_degrees())

        for key in self._jdd:
            self._jdd[key] = self._evaluate_prob_of_joint_degree(ast.literal_eval(key))
        
        self.normalise_jdd()

    def _draw_from_analytical_joint(self)->np.ndarray:
        '''Draw `self._n_samples' from a marginal distribution `self._arr_fp' along each dimension.
        :returns jds: samples from marginal distributions'''
        ret = []
        for i in range(len(self._hi_lo_degree_bounds)):
            kmin,kmax = self._hi_lo_degree_bounds[i]
            ks = [k for k in range(kmin, kmax+1)]                   # possible degrees
            pks = [self._arr_fp[i](k) for k in ks]                  # degree weights 
            ret.append(random.choices(ks,pks,k=self._n_samples))    # sample this dimension
        return np.column_stack(ret)                                 # return sampled degrees

    def _create_jdd_by_sampling(self)->None:
        '''Draws jds samples from analytical marginal functions and converts to 
        a joint degree distribution. '''
        self.convert_jds_to_jdd(self._draw_from_analytical_joint())

class JDD_split_K_model(JDD_Interface):
    '''An overall degree distribution is split with probabilities of creating each motif.'''

    def __init__(self, fp : Callable, 
                       motif_sizes : List[int],
                       probs : List[float],
                       kmin : int, 
                       kmax : int)->None:
        
        super().__init__(motif_sizes)
        self._fp = fp
        self._probs = probs
        self._kmin = kmin
        self._kmax = kmax
        
        self.create_jdd()
        
    def create_jdd(self)->None:
        '''Creates self._jdd using the split degree model'''
        for k in range(self._kmin,self._kmax):
            self.resolve_degree(k, self._fp(k))
        self.normalise_jdd()

    def get_valid_joint_degrees(self, remaining_degree : int, topology : int)->List[_JOINT_DEGREE]:
        '''Returns a list of tuples by recursion. Only an ordered list of cliques are currently supported.
        
        :param remaining_degree: current free edges that can be partitioned
        :param topology: column index of joint degree tuple
        
        :returns list of joint degrees'''
        if topology  == 1:
            yield [remaining_degree]
        else:
            for i in range( 0, remaining_degree//topology+1 ):
                for row in self.get_valid_joint_degrees(remaining_degree-i*topology, topology-1):
                    yield row+[i]
    
    def calc_prob_of_joint_degree(self, jd : _JOINT_DEGREE)->float:
        '''calculates the probability of a joint degree from the input params.
        
        :param jd: joint degree
        :returns probability: float value'''

        prod : float = 1.0
        for i, degree in enumerate(jd):
            prod *= pow(self._probs[i], (i+1)*degree)
        return prod

    def resolve_degree(self, k : int, prob_overall_k : float)->None:
        '''creates all valid joint degrees for overall k and their probability and updates self._jdd.
        
        :param k: overall degree
        :param prob_overall_k: float value'''
        
        # get a list of valid joint degrees
        valid_tuples : List[_JOINT_DEGREE] = list(self.get_valid_joint_degrees(k,len(self._probs)))

        # calculate the probability of each joint degree tuple
        probabilities = []
        for jd in valid_tuples:
            probabilities.append(self.calc_prob_of_joint_degree(jd))

        # normalise the probabilities to unity
        total = sum(probabilities)

        for i in range(len(probabilities)):
            probabilities[i] /= total

        # add each tuple to self._jdd weighted by probability of overall degree k
        for i, jd in enumerate(valid_tuples):
            self._jdd[str(jd)] = prob_overall_k * probabilities[i]

class JDD_delta_model(JDD_split_K_model):
    '''A distribution of single edges apart from a specified target degree. For instance
    a distribution of degrees in a single dimension (2-cliques) with zeros for all other motif
    type counts, apart from when k=target. '''
    def __init__(self, fp : Callable, 
                       motif_sizes : List[int],
                       probs : List[float],
                       target_k : int, 
                       kmin : int, 
                       kmax : int)->None:

        self._target_k = target_k
        super().__init__(fp, motif_sizes, probs, kmin, kmax)
        
    def create_jdd(self)->None:
        '''Creates self._jdd using the degree delta model. '''
        for k in range(self._kmin,self._kmax):
            zeros = [0]*len(self._motif_sizes)
            if k != self._target_k:
                zeros[0] = k
                self._jdd[str(tuple(zeros))] = self._fp(k)
            else:
                self.resolve_degree(k, self._fp(k))

        self.normalise_jdd()

class JDD_clique_cover(JDD_Interface):
    '''Creates self._jdd from a list of cliques in the network. The sizes of the 
    cliques can be obtained from the self._motif_sizes member.'''

    def __init__(self, C : _COVER):
        ''':param C: clique cover'''
        self._cover = C
        super().__init__(sorted(list(set([len(c) for c in C]))))
        self.create_jdd()

    def create_jdd(self)->None:

        node_ids = list(set([node for clique in self._cover for node in clique]))

        zero_index = 0
        if min(node_ids) != zero_index:
            zero_index = 1

        largest_clique = len(max(self._cover, key = len))

        jds = []
        for _ in range(len(node_ids)):
            jd = [0] * largest_clique
            jds.append(jd)

        for c in self._cover:

            clique_size = len(c)
            for node in c:
                jds[node - zero_index][clique_size-1] += 1

        # iterate each column of the jds and record the index if all zeros
        indxs = [i for i, top in enumerate(zip(*jds)) if not any(top)]
        
        # use the indexes of the zero columns to remove
        for i in indxs:
            for jd in jds:
                del jd[i]

        # convert jds to jdd
        self.convert_jds_to_jdd(jds)
        


