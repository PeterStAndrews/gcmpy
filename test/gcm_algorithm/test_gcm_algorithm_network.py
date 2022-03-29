

import unittest

from gcmpy.motif_generators.clique_motif import clique_motif
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import JointDegreeManual
from gcmpy.names.joint_degree_names import JointDegreeNames
from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames
from gcmpy.gcm_algorithm.gcm_algorithm_network import GCMAlgorithmNetwork

NETWORK_SIZE: int = 100000

class GCMAlgorithmNetworkTest(unittest.TestCase):

    def test_single_topology(self):

        params = {}
        params[JointDegreeNames.JDD] = {(1,) : 0.2, (2,) : 0.5, (3,) : 0.1, (5,) : 0.2}
        params[JointDegreeNames.MOTIF_SIZES] = [2]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices : int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params[GCMAlgorithmNames.MOTIF_SIZES] = [2]
        params[GCMAlgorithmNames.EDGE_NAMES] = ['2-clique']
        params[GCMAlgorithmNames.BUILD_FUNCTIONS] = [clique_motif]
        g = GCMAlgorithmNetwork(
                params
            ).random_clustered_graph(jds)


        num_expected_edges = sum([k[0] for k in jds]) / 2.0

        # assert that random graph has +- 0.5% of expected edges
        tolerance = num_expected_edges / 1000.0

        self.assertTrue(g._G.order() == NETWORK_SIZE)
        self.assertTrue(num_expected_edges - tolerance <= len(g._G.edges) <= num_expected_edges + tolerance)

    def test_two_topologies(self):

        params = {}
        params[JointDegreeNames.JDD] = {(1,0) : 0.2, (2,1) : 0.5, (3,0) : 0.1, (5,1) : 0.2}
        params[JointDegreeNames.MOTIF_SIZES] = [2,3]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices : int = NETWORK_SIZE 
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params[GCMAlgorithmNames.MOTIF_SIZES] = [2,3]
        params[GCMAlgorithmNames.EDGE_NAMES] = ['2-clique','3-clique']
        params[GCMAlgorithmNames.BUILD_FUNCTIONS] = [clique_motif,clique_motif]
        g = GCMAlgorithmNetwork(
                params
            ).random_clustered_graph(jds)

        
        num_2_clique_edges = sum([k[0] for k in jds]) / 2.0
        num_3_clique_edges = sum([k[1] for k in jds])
        num_expected_edges = num_2_clique_edges + num_3_clique_edges
        
        # assert that random graph has +- 0.5% of expected edges
        tolerance = num_expected_edges / 100.0
        
        self.assertTrue(g._G.order() == NETWORK_SIZE)
        num_edges: int = len(g._G.edges)
        self.assertTrue(
            num_expected_edges-tolerance <= num_edges <= num_expected_edges+tolerance
        )