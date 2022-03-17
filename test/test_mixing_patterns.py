

import unittest
from gcmpy.tools.joint_excess_joint_degree import MotifMixingPatterns
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import JointDegreeManual
from gcmpy.motif_generators.clique_motif import clique_motif
from gcmpy.gcm_algorithm.gcm_algorithm_network import GCMAlgorithmNetwork

NETWORK_SIZE: int = 100000

class MotifMixingPatternsTest(unittest.TestCase):

    def test_two_topologies(self):

        params = {}
        params["jdd"] = {(1,0) : 0.2, (2,1) : 0.5, (3,0) : 0.1, (5,1) : 0.2}
        params["motif_sizes"] = [2,3]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices : int = NETWORK_SIZE 
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params["motif_sizes"] = [2,3]
        params["edge_names"] = ['2-clique','2-clique']
        params["build_functions"] = [clique_motif,clique_motif]
        g = GCMAlgorithmNetwork(
                params
            ).random_clustered_graph(jds)

        params = {}
        params['network'] = g._G
        params['jdd'] = {(5,1) : 1/3, (3,2) : 1/3, (1,3) : 1/3}
        params['edge_names'] = ['2-clique','2-clique']
        C = MotifMixingPatterns(params)

        ejks_1 = C.get_ejks()

        # call again without supplying jdd
        params = {}
        params['network'] = g._G
        params['edge_names'] = ['2-clique','2-clique']
        C = MotifMixingPatterns(params)

        ejks_2 = C.get_ejks()

        self.assertTrue(ejks_1==ejks_2)

    def test_three_topologies(self):

        params = {}
        params["jdd"] = {(1,0,0) : 0.2, (2,1,1) : 0.5, (3,0,1) : 0.1, (5,1,0) : 0.2}
        params["motif_sizes"] = [2,3,2]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices : int = NETWORK_SIZE 
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params["edge_names"] = ['2-clique-blue', '3-clique','2-clique-red']
        params["build_functions"] = [clique_motif,clique_motif,clique_motif]
        params["motif_sizes"] = [2,3,2]
        g = GCMAlgorithmNetwork(
                params
            ).random_clustered_graph(jds)

        params = {}
        params["network"] = g._G
        params["edge_names"] = ['2-clique-blue', '3-clique','2-clique-red']
        C = MotifMixingPatterns(params)

        self.assertTrue(len(C.get_ejks())==3)