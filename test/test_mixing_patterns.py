

import unittest
from gcmpy.tools.joint_excess_joint_degree import MotifMixingPatterns
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import JointDegreeManual
from gcmpy.motif_generators.clique_motif import clique_motif
from gcmpy.gcm_algorithm.gcm_algorithm import GCMAlgorithm

NETWORK_SIZE: int = 100000

class MotifMixingPatternsTest(unittest.TestCase):

    def test_two_topologies(self):

        params = {}
        params["jdd"] = {(1,0) : 0.2, (2,1) : 0.5, (3,0) : 0.1, (5,1) : 0.2}
        params["motif_sizes"] = [2,3]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices : int = NETWORK_SIZE 
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        edge_names: list[str] = ['2-clique', '3-clique']
        build_functions: list[callable] = [clique_motif,clique_motif]
        g = GCMAlgorithm(
                params["motif_sizes"], build_functions, edge_names
            ).random_clustered_graph(jds)

        params = {}
        params['network'] = g._G
        params['jdd'] = {(5,1) : 1/3, (3,2) : 1/3, (1,3) : 1/3}
        params['topology_names'] = edge_names
        C = MotifMixingPatterns(params)

        ejks_1 = C.get_ejks()

        # call again without supplying jdd
        params = {}
        params['network'] = g._G
        params['topology_names'] = edge_names
        C = MotifMixingPatterns(params)

        ejks_2 = C.get_ejks()

        self.assertTrue(ejks_1==ejks_2)
