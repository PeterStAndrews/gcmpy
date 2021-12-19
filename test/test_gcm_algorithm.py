

from gcmpy import *
import unittest
import time

class GCM_algorithm_Test(unittest.TestCase):

    def test_manual_single_topology(self):
        
        # valid input data for manual entry
        motif_sizes = [2]
        jdd = {'(1,)' : 0.2, '(2,)' : 0.5, '(3,)' : 0.1, '(5,)' : 0.2}

        # create degree distribution
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        # draw a sample
        n_vertices : int = NETWORK_SIZE
        jds = DegreeDistObj.sample_JDS(n_vertices)

        # create graph 
        build_functions = [clique_motif]
        g = GCM_algorithm(motif_sizes, build_functions).random_clustered_graph(jds)
        
    def test_marginal_JDD_single_topology(self):

        motif_sizes = [2]                            # 2-cliques
        powerlaw_exponent = 2.5                      # powerlaw_exponent
        fp_array = [power_law(powerlaw_exponent)]    # array of marginals
        kmax = 150000                                # largest degree
        kmin = 1                                     # smallest degree

        # create joint degree distribution object
        DegreeDistObj = JDD_marginals(fp_array, motif_sizes, [(kmin,kmax)])

        # draw a sample
        n_vertices : int = NETWORK_SIZE 
        jds = DegreeDistObj.sample_JDS(n_vertices)

        # create graph 
        build_functions = [clique_motif]
        g = GCM_algorithm(motif_sizes, build_functions).random_clustered_graph(jds)

    def test_split_k_two_topologies(self):

        motif_sizes = [2,3]                    
        probs = [0.8,0.2]   # probability that an edge is a 2- or 3-clique

        powerlaw_exponent = 2.5                     
        fp = power_law(powerlaw_exponent)   # overall degree distribution

        kmax = 1000                                  
        kmin = 1                                     

        DegreeDistObj = JDD_split_K_model(fp, motif_sizes, probs, kmin, kmax)

        n_vertices : int = NETWORK_SIZE
        jds = DegreeDistObj.sample_JDS(n_vertices)
        
        # create graph 
        build_functions = [clique_motif, clique_motif]
        g = GCM_algorithm(motif_sizes, build_functions).random_clustered_graph(jds)

class ResampleJDS_Test(unittest.TestCase):

    def test_manual_single_topology(self):

         # valid input data for manual entry
        motif_sizes = [2]
        jdd = {'(1,)' : 0.2, '(2,)' : 0.5, '(3,)' : 0.1, '(5,)' : 0.2}

        # create degree distribution
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        # draw a sample
        n_vertices : int = NETWORK_SIZE
        jds = DegreeDistObj.sample_JDS(n_vertices)

        # create graph 
        num_networks = 10
        build_functions = [clique_motif]
        g = ResampleJDS(num_networks, motif_sizes, build_functions).random_clustered_graph(jds)
        