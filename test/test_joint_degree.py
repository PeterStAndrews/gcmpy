
# Test the behaviour of the joint degree distribution 
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

from gcmpy import *
import unittest
import time

class JDD_manual_Test(unittest.TestCase):

    def test_manual_JDD_single_topology(self):
        
        # valid input data for manual entry
        motif_sizes = [2]
        jdd = {'(1,)' : 0.2, '(2,)' : 0.5, '(3,)' : 0.1, '(5,)' : 0.2}

        # create degree distribution
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)
        
    def test_manual_JDD_single_topology_performance(self):

        # start timer
        t0 = time.time()

        # valid input data for manual entry
        motif_sizes = [2]
        jdd = {'(1,)' : 0.2, '(2,)' : 0.5, '(3,)' : 0.1, '(5,)' : 0.2}

        # create degree distribution
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)
        t1 = time.time()

        total = t1-t0
        self.assertGreaterEqual(18, total)

    def test_manual_JDD_two_topologies(self):
        
        # valid input data for manual entry
        motif_sizes = [2,3]
        jdd = {'(1,0)' : 0.2, '(2,1)' : 0.5, '(3,0)' : 0.1, '(5,1)' : 0.2}

        # create degree distribution
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

class JDD_empirical_data_Test(unittest.TestCase):

    def test_empirical_JDD_single_topology_list(self):

        random_degrees = np.random.randint(1,20,5000)
        jds = [ [k] for k in random_degrees]

        motif_sizes = [2]                    
        DegreeDistObj = JDD_empirical_data(jds, motif_sizes)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

    def test_empirical_JDD_single_topology_tuple(self):

        jds = [(k,) for k in np.random.randint(1,20,5)]
        motif_sizes = [2]
        DegreeDistObj = JDD_empirical_data(jds, motif_sizes)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)
        
    def test_empirical_JDD_two_topology_list(self):

        jds = [[k,k+1] for k in np.random.randint(1,20,5)]
        motif_sizes = [2,3]
        DegreeDistObj = JDD_empirical_data(jds, motif_sizes)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

    def test_empirical_JDD_two_topology_tuple(self):

        jds = [(k,k+1) for k in np.random.randint(1,20,5)]
        motif_sizes = [2,3]
        DegreeDistObj = JDD_empirical_data(jds, motif_sizes)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

class JDD_marginal_Test(unittest.TestCase):

    def test_marginal_JDD_single_topology(self):

        motif_sizes = [2]                    # 2-cliques
        mean_degree = 2.5                    # mean poisson degree
        fp_array = [poisson(mean_degree)]    # array of marginals
        kmax = 15                            # largest degree
        kmin = 0                             # smallest degree

        # create joint degree distribution object
        DegreeDistObj = JDD_marginals(fp_array, motif_sizes, [(kmin,kmax)])
        
        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

    def test_marginal_JDD_two_topologies(self):

        # allow 2 & 3-cliques
        motif_sizes = [2,3]
        
        # mean poisson degrees
        mean_degree_2_clique = 2.5
        mean_degree_3_clique = 0.5
        
        # array of marginal distribution functions
        fp_array = [poisson(mean_degree_2_clique),
                    poisson(mean_degree_3_clique)]

        # kmin and kmax in each dimension
        kmax_2_clique = 15                            
        kmin_2_clique = 0      

        kmax_3_clique = 5                           
        kmin_3_clique = 0   

        kmin_max = [(kmin_2_clique,kmax_2_clique),(kmin_3_clique,kmax_3_clique)]                           

        # create joint degree distribution object
        DegreeDistObj = JDD_marginals(fp_array, motif_sizes, kmin_max)
        
        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

    def test_marginal_JDD_two_topologies_sample(self):

        # allow 2 & 3-cliques
        motif_sizes = [2,3]
        
        # mean poisson degrees
        mean_degree_2_clique = 2.5
        mean_degree_3_clique = 0.5
        
        # array of marginal distribution functions
        fp_array = [poisson(mean_degree_2_clique),
                    poisson(mean_degree_3_clique)]

        # kmin and kmax in each dimension
        kmax_2_clique = 15                            
        kmin_2_clique = 0      

        kmax_3_clique = 5                           
        kmin_3_clique = 0   

        kmin_max = [(kmin_2_clique,kmax_2_clique),(kmin_3_clique,kmax_3_clique)]                           

        # create joint degree distribution object by sampling
        DegreeDistObj = JDD_marginals(fp_array, motif_sizes, kmin_max,
                                      use_sampling=True,n_samples=int(1e5))

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

    def test_marginal_JDD_many_topologies_sample(self):

        # allow up to n-cliques
        n = 10
        motif_sizes = list(range(2,n+1))
        
        # mean poisson degrees (all same for ease)
        mean_degree_clique = 2
        
        # array of marginal distribution functions (all same for ease)
        fp_array = [poisson(mean_degree_clique)] * (n+1-2)

        # kmin and kmax (all same for ease)
        kmax = 10                            
        kmin= 0   

        kmin_max = [(kmin,kmax)] * (n+1-2)                        

        # create joint degree distribution object by sampling
        DegreeDistObj = JDD_marginals(fp_array, motif_sizes, kmin_max,
                                      use_sampling=True,n_samples=int(1e5))

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

    def test_marginal_JDD_powerlaw_topology(self):

        motif_sizes = [2]                            # 2-cliques
        powerlaw_exponent = 2.5                      # powerlaw_exponent
        fp_array = [power_law(powerlaw_exponent)]    # array of marginals
        kmax = 150000                                # largest degree
        kmin = 1                                     # smallest degree

        # create joint degree distribution object
        DegreeDistObj = JDD_marginals(fp_array, motif_sizes, [(kmin,kmax)])

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

    def test_marginal_JDD_scale_free_cut_off_topology(self):

        motif_sizes = [2]                            # 2-cliques
        powerlaw_exponent = 2.5                      # powerlaw_exponent
        degree_cut_off = 25                          # degree cut off 

                                                     # array of marginals            
        fp_array = [scale_free_cut_off(powerlaw_exponent,degree_cut_off)]    
        kmax = 1500                                  # largest degree
        kmin = 1                                     # smallest degree

        # create joint degree distribution object
        DegreeDistObj = JDD_marginals(fp_array, motif_sizes, [(kmin,kmax)])

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

class JDD_split_K_Test(unittest.TestCase):

    def test_split_k_two_topologies(self):

        motif_sizes = [2,3]                    
        probs = [0.8,0.2]   # probability that an edge is a 2- or 3-clique

        powerlaw_exponent = 2.5                     
        fp = power_law(powerlaw_exponent)   # overall degree distribution

        kmax = 1000                                  
        kmin = 1                                     

        DegreeDistObj = JDD_split_K_model(fp, motif_sizes, probs, kmin, kmax)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

class JDD_delta_Test(unittest.TestCase):

    def test_delta_two_topologies(self):
        
        motif_sizes = [2,3]                    
        probs = [0.8,0.2]   # probability that an edge is a 2- or 3-clique

        powerlaw_exponent = 2.5                     
        fp = power_law(powerlaw_exponent)   # overall degree distribution

        kmax = 1000                                  
        kmin = 1                                     

        target_k = 3
        DegreeDistObj = JDD_delta_model(fp, motif_sizes, probs, target_k, kmin, kmax)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)

def EECC_Test_util_network()->EECC:

    G = EECC()

    G.add_edge((1, 2))
    G.add_edge((1, 14))
    G.add_edge((2, 4))
    G.add_edge((2, 13))
    G.add_edge((2, 14))
    G.add_edge((3, 4))
    G.add_edge((3, 5))
    G.add_edge((4, 5))
    G.add_edge((4, 13))
    G.add_edge((4, 14))
    G.add_edge((6, 7))
    G.add_edge((6, 13))
    G.add_edge((7, 8))
    G.add_edge((7, 13))
    G.add_edge((8, 9))
    G.add_edge((8, 13))
    G.add_edge((9, 10))
    G.add_edge((9, 11))
    G.add_edge((9, 13))
    G.add_edge((10, 11))
    G.add_edge((11, 12))
    G.add_edge((12, 13))
    G.add_edge((13, 14))

    return G

class JDD_cover_Test(unittest.TestCase):

    def test_EECC_cover(self):
        
        G = EECC_Test_util_network()

        G.set_max_clique_size(4)
        EECC_cover = G.get_EECC()

        DegreeDistObj = JDD_clique_cover(EECC_cover)

        n_vertices : int = NETWORK_SIZE
        DegreeDistObj.sample_JDS(n_vertices)









if __name__ == '__main__':
    unittest.main()