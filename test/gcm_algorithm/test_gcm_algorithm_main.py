

import unittest

from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import JointDegreeManual
from gcmpy.motif_generators.clique_motif import clique_motif
from gcmpy.gcm_algorithm.gcm_algorithm_main import GeneralisedConfigurationModel
from gcmpy.gcm_algorithm.gcm_algorithm_types import GCMAlgorithmTypes
from gcmpy.network.network import Network
from gcmpy.network.edge_list import LightWeightEdgeList

NETWORK_SIZE = 100000

class GeneralisedConfigurationModelTest(unittest.TestCase):

    def test_main_entrypoint(self):

        params = {}
        params["jdd"] = {(1,) : 0.2, (2,) : 0.5, (3,) : 0.1, (5,) : 0.2}
        params["motif_sizes"] = [2]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices: int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params["GCM_type"] = GCMAlgorithmTypes.FAST
        params["motif_sizes"] = [2]
        params["edge_names"] = ['2-clique']
        params["build_functions"] = [clique_motif]

        # test edge list
        g: LightWeightEdgeList = GeneralisedConfigurationModel.load_gcm_algorithm(
            params
        ).random_clustered_graph(jds)

        # test network
        params["GCM_type"] = GCMAlgorithmTypes.NETWORK
        g: Network = GeneralisedConfigurationModel.load_gcm_algorithm(
            params
        ).random_clustered_graph(jds)
