
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

        n_vertices : int = 10000
        DegreeDistObj.sample_JDS(n_vertices)
        
    def test_manual_JDD_single_topology_performance(self):

        # start timer
        t0 = time.time()

        # valid input data for manual entry
        motif_sizes = [2]
        jdd = {'(1,)' : 0.2, '(2,)' : 0.5, '(3,)' : 0.1, '(5,)' : 0.2}

        # create degree distribution
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        n_vertices : int = 10000
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

        n_vertices : int = 1000
        DegreeDistObj.sample_JDS(n_vertices)

class JDD_empirical_data_Test(unittest.TestCase):

    def test_empirical_JDD_single_topology_list(self):

        random_degrees = np.random.randint(1,20,5000)
        jds = [ [k] for k in random_degrees]

        motif_sizes = [2]                    
        DegreeDistObj = JDD_empirical_data(jds, motif_sizes)

        n_vertices : int = 10000
        DegreeDistObj.sample_JDS(n_vertices)

    def test_empirical_JDD_single_topology_tuple(self):

        jds = [(k,) for k in np.random.randint(1,20,5)]
        motif_sizes = [2]
        DegreeDistObj = JDD_empirical_data(jds, motif_sizes)

        n_vertices : int = 10000
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
        
        n_vertices : int = 10000
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
        
        n_vertices : int = 10000
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

        n_vertices : int = 10000
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

        n_vertices : int = 10000
        DegreeDistObj.sample_JDS(n_vertices)

    def test_marginal_JDD_powerlaw_topology(self):

        motif_sizes = [2]                            # 2-cliques
        powerlaw_exponent = 2.5                      # powerlaw_exponent
        fp_array = [power_law(powerlaw_exponent)]    # array of marginals
        kmax = 150000                                # largest degree
        kmin = 1                                     # smallest degree

        # create joint degree distribution object
        DegreeDistObj = JDD_marginals(fp_array, motif_sizes, [(kmin,kmax)])

        n_vertices : int = 10000
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

        n_vertices : int = 10000
        DegreeDistObj.sample_JDS(n_vertices)



if __name__ == '__main__':
    unittest.main()