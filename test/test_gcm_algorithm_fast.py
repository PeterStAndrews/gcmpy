
import unittest


from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import JointDegreeManual
from gcmpy.motif_generators.clique_motif import clique_motif
from gcmpy.gcm_algorithm.gcm_algorithm_fast import FastGCMAlgorithm

NETWORK_SIZE: int = 100000

class GCMAlgorithmFastTest(unittest.TestCase):

    def test_single_topology(self):

        params = {}
        params["jdd"] = {(1,) : 0.2, (2,) : 0.5, (3,) : 0.1, (5,) : 0.2}
        params["motif_sizes"] = [2]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices: int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params["motif_sizes"] = [2]
        params["edge_names"] = ['2-clique']
        params["build_functions"] = [clique_motif]
        g = FastGCMAlgorithm(
                params
            ).random_clustered_graph(jds)

        num_edges = sum([k[0] for k in jds]) / 2.0

        self.assertTrue(num_edges - 2 <= len(g._edge_list) <= num_edges + 2)

    def test_two_topologies(self):

        params = {}
        params["jdd"] = {(1,0) : 0.2, (2,1) : 0.5, (3,0) : 0.1, (5,1) : 0.2}
        params["motif_sizes"] = [2,3]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices : int = NETWORK_SIZE 
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params["motif_sizes"] = [2,3]
        params["edge_names"] = ['2-clique','3-clique']
        params["build_functions"] = [clique_motif,clique_motif]
        g = FastGCMAlgorithm(
                params
            ).random_clustered_graph(jds)

        tolerance = 20
        num_2_clique_edges = sum([k[0] for k in jds]) / 2.0
        num_3_clique_edges = sum([k[1] for k in jds])
        num_expected_edges = num_2_clique_edges + num_3_clique_edges
        
        x: int = len(g._edge_list)
        self.assertTrue(
            num_expected_edges-tolerance <= x <= num_expected_edges+tolerance
        )