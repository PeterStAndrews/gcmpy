

from gcmpy import *
import unittest
import cProfile
import functools
import pstats

def profile(func):

    @functools.wraps(func)
    def inner(*args, **kwargs):
        profiler = cProfile.Profile()
        profiler.enable()
        try:
            retval = func(*args, **kwargs)
        finally:
            profiler.disable()
            with open('profile.out', 'w') as profile_file:
                stats = pstats.Stats(profiler, stream=profile_file)
                stats.print_stats()
        return retval

    return inner

@profile
def _profle_gcm(gcm : GCM_Network_Generator):
    gcm.create_gcm_networks()

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
        
    def test_node_edge_attributes(self):
        
        # valid input data for manual entry
        motif_sizes = [2]
        jdd = {'(1,)' : 0.2, '(2,)' : 0.5, '(3,)' : 0.1, '(5,)' : 0.2}

        edge_topology_names = ['2-clique']

        # create degree distribution
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        # draw a sample
        n_vertices : int = NETWORK_SIZE
        jds = DegreeDistObj.sample_JDS(n_vertices)

        # create graph 
        build_functions = [clique_motif]
        model = GCM_algorithm(motif_sizes,
                              build_functions,
                              edge_topology_names).random_clustered_graph(jds)


    def test_marginal_JDD_single_topology(self):

        motif_sizes = [2]                            # 2-cliques
        powerlaw_exponent = 2.5                      # powerlaw_exponent
        fp_array = [power_law(powerlaw_exponent)]    # array of marginals
        kmax = 1500                                  # largest degree
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

class GCM_Network_Generator_Test(unittest.TestCase):

    def test_manual_single_topology(self):

         # create the degree distribution object
        motif_sizes = [2]
        jdd = {'(1,)' : 0.2, '(2,)' : 0.5, '(3,)' : 0.1, '(5,)' : 0.2}
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        num_networks = 10                   # number of networks to create
        n_vertices : int = NETWORK_SIZE     # number of vertices in the networks
        build_functions = [clique_motif]    # motif build functions
        resample_jds_every : int = 1000     # resample a new JDS every n steps

        # create GCM object with parameters
        gcm = GCM_Network_Generator(num_networks,
                                    motif_sizes,
                                    build_functions,
                                    n_vertices,
                                    resample_jds_every)

        # insert degree distribution object
        gcm.insert_jdd_generator(DegreeDistObj)

        # create the networks and return the results
        gcm.create_gcm_networks()

    def test_profile(self):
        
        motif_sizes = [2]
        jdd = {'(1,)' : 0.2, '(2,)' : 0.5, '(3,)' : 0.1, '(5,)' : 0.2}
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        num_networks = 1                    # number of networks to create
        n_vertices : int = NETWORK_SIZE     # number of vertices in the networks
        build_functions = [clique_motif]    # motif build functions
        resample_jds_every : int = 1        # resample a new JDS every n steps

        # create GCM object with parameters
        gcm = GCM_Network_Generator(num_networks,
                                    motif_sizes,
                                    build_functions,
                                    n_vertices,
                                    resample_jds_every)

        # insert degree distribution object
        gcm.insert_jdd_generator(DegreeDistObj)

        # create the networks and profile the code.
        _profle_gcm(gcm)
       

    def test_generator(self):

        # parameters for network creation (gcmpy)
        motif_sizes = [2]                            # 2-cliques

        NETWORK_SIZE = 6000
        build_functions = [clique_motif]

        n_repeats = 2
        motif_sizes = [2]
        jdd = {'(1,)' : 0.2, '(2,)' : 0.5, '(3,)' : 0.1, '(5,)' : 0.2}
        DegreeDistObj = JDD_manual(jdd, motif_sizes)

        n_sample_mod = 1
        NetGenerator = GCM_Network_Generator(n_repeats, 
                                                   motif_sizes, 
                                                   build_functions, 
                                                   NETWORK_SIZE, 
                                                   n_sample_mod)
        
        # insert joint degree distribution object
        NetGenerator.insert_jdd_generator(DegreeDistObj)
        
        # create list of networks
        Gs = NetGenerator.create_gcm_networks()
        