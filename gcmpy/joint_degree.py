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

import numpy as np
import random
import ast
from itertools import product
from typing import Tuple, Dict, List, Callable, TypeVar
from collections import Counter

_JDD   = TypeVar('_JDD', bound=Dict[str,float])  
_JDS   = TypeVar('_JDS', bound=List[List[int]])          
_JOINT_DEGREE = TypeVar('_JOINT_DEGREE', bound=Tuple[int,...])

class JDD_Interface(object):
    '''Joint degree distribution interface. Subclasses will define how self._jdd
    is created, and all implementations should do so upon construction.'''

    def __init__( self, modulo ):
        self._jdd : _JDD
        self._modulo : List[int] = modulo

    def _sample_JDS( self, N : int )->_JDS:
        '''Samples N joint motif sequence from self._jdd object. Intended to 
        be a private method, only to be called by sample_JDS. '''
        keys = np.array(list(self._jdd.keys()))
        values = np.array(list(self._jdd.values()))
        choice_list_indices = np.random.choice(len(keys), N, replace=True, p=values)
        return keys[choice_list_indices] 

    def sample_JDS( self, N : int )->None:  
        '''Method to ensure the sum of the motifs in each dimension mod modulo[i] (which is the 
        number of vertices in each motif) is zero by adding motifs to the motif list. 
        This ensures the JDS is graphic.'''
        # raw joint degree sequence
        jds = list(map(ast.literal_eval,self._sample_JDS(N)))
        
        # ensure that the sum of the numbers is divisible by the motif size
        ntops = list(map(sum, zip(*jds)))
        for i, ntop in enumerate(ntops):
            if ntop % self._modulo[i] != 0:
                # if not, round up to add one more motif to the network
                for j in range(self._modulo[i] - ntop % self._modulo[i]):
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
        '''Convert a joint degree sequence to self._jdd .'''
        n_samples = len(jds)
        self._jdd = dict((k,v/n_samples) for k, v in Counter(str(tuple(e)) for e in list(jds)).items())
 
class JDD_manual(JDD_Interface):
    '''Specify self._jdd by hand.''' 
    def __init__(self, jdd : _JDD, modulo : List[int])->None:
        self._jdd = jdd
        super().__init__(modulo)

class JDD_empirical_data(JDD_Interface):
    '''An empirical joint degree sequence is used to create self._jdd '''

    def __init__(self, jds : _JDS, modulo : List[int])->None:
        self.convert_jds_to_jdd(jds)
        super().__init__(modulo)

class JDD_joint_function(JDD_Interface):
    '''Multivariate function to evaluate the probability of given joint degree from an analytical source.
    Note, callable self._fp must accept a joint degree tuple and return a float.'''
    def __init__(self, fp : Callable,
                       modulo : List[int],
                       hi_lo_degree_bounds : Tuple[int,int],
                       use_sampling : bool = False,
                       n_samples : int = 1e5
                       ):
        self._fp = fp
        self._hi_lo_degree_bounds = hi_lo_degree_bounds  
        super().__init__(modulo)
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
    by setting `use_sampling' which draws `n_samples' weighted samples from each marginal function.'''
    def __init__(self, arr_fp : List[Callable],
                       modulo : List[int],
                       hi_lo_degree_bounds : List[Tuple[int,int]],
                       use_sampling : bool = False,
                       n_samples : int = 1e5
                       ):
        self._arr_fp = arr_fp
        self._hi_lo_degree_bounds = hi_lo_degree_bounds

        if not use_sampling:
            self._create_jdd_directly()
        else:
            self._n_samples = n_samples
            self._create_jdd_by_sampling()

        super().__init__(modulo)

    def _generate_all_joint_degrees(self)->List[_JOINT_DEGREE] :
        '''Generate all possible joint degrees from a range of min/max degree of each motif.'''
        ks = []
        for kmin, kmax in self._hi_lo_degree_bounds:
            ks.append([k for k in range(kmin,kmax)])
        return list(product(*ks))
    
    def _evaluate_prob_of_joint_degree(self, joint_degree : List[_JOINT_DEGREE])->float:
        '''Evaluate the joint probability of a joint degree using function callbacks'''
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
        
        self.normalise_jdd(self)

    def _draw_from_analytical_joint(self)->np.ndarray:
        '''Draw `self._n_samples' from a marginal distribution `self._arr_fp' along each dimension.'''
        ret = []
        for i in range(len(self._hi_lo_degree_bounds)):
            kmin,kmax = self._hi_lo_degree_bounds[i]
            ks = [k for k in range(kmin, kmax+1)]                   # possible degrees
            pks = [self._arr_fp[i](k) for k in ks]                  # degree weights 
            ret.append(random.choices(ks,pks,k=self._n_samples))    # sample this dimension
        return np.column_stack(ret)                                 # return sampled degrees

    def _create_jdd_by_sampling(self)->None:
        self.convert_jds_to_jdd(self._draw_from_analytical_joint())

class JDD_split_K_model(JDD_Interface):
    '''An overall degree distribution is split with probabilities of creating each motif.'''

    def __init__(self, fp : Callable, 
                       modulo : List[int],
                       probs : List[float],
                       num_edges : List[int], 
                       kmin : int, 
                       kmax : int)->None:

        self._fp = fp
        self._probs = probs
        self._hi_lo = (kmin,kmax)
        self._num_edges = num_edges                 # number of edges per vertex per motif
        self.create_jdd()
        super().__init__(modulo)

    def create_jdd(self)->None:
        '''Creates self._jdd using the split degree model'''
        for k in range(self._kmin,self._kmax):
            self.resolve_degree(k, self._fp(k))

    def get_valid_joint_degrees(self, target : int, column : int)->List[_JOINT_DEGREE]:
        '''Returns a list of tuples by recursion. Only an ordered list of cliques are currently supported.'''
        if column == 1:
            yield [target]
        else:
            for i in range( 0, target//column+1 ):
                for row in self.get_valid_joint_degrees( target-i*column, column-1 ):
                    yield row+[i]
    
    def calc_prob_of_joint_degree(self, jd : _JOINT_DEGREE)->float:
        '''calculates the probability of a joint degree from the input params.'''
        prod : float = 1.0
        for i, degree in enumerate(jd):
            prod *= pow(self._probs[i],self._num_edges[i]*degree)
        return prod

    def resolve_degree(self, k : int, prob_overall_k : float)->None:
        '''creates all valid joint degrees for overall k and their probability and updates self._jdd'''
        
        # get a list of valid joint degrees
        valid_tuples : List[_JOINT_DEGREE] = list(self.get_valid_joint_degrees(k))

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
                       modulo : List[int],
                       probs : List[float],
                       num_edges : List[int],
                       target_k : int, 
                       kmin : int, 
                       kmax : int)->None:

        self._target_k = target_k
        super().__init__(fp, modulo, probs, num_edges, kmin, kmax)

    def create_jdd(self)->None:
        '''Creates self._jdd using the degree delta model. '''
        for k in range(self._kmin,self._kmax):
            zeros = [0]*len(self._modulo)
            if k != self._target_k:
                zeros[0] = k
                self._jdd[str(tuple(zeros))] = self._fp(k)
            else:
                self.resolve_degree(k, self._fp(k))
