

import unittest

from gcmpy.tools.joint_degree_distribution_from_network import JointDegreeDistributionFromNetwork
from gcmpy.joint_degree.joint_degree_distribution import JointDegreeDistribution
from gcmpy.joint_degree.joint_degree_types import JointDegreeType
from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames
from gcmpy.names.joint_degree_names import JointDegreeNames

from gcmpy.motif_generators.clique_motif import clique_motif
from gcmpy.gcm_algorithm.gcm_algorithm_network import GCMAlgorithmNetwork

NETWORK_SIZE: int = 100000

class JointDegreeDistributionFromNetworkTest(unittest.TestCase):

    def test_jdd_from_net(self):

        jdd = {(1,2) : 0.2, (2,0) : 0.5, (3,1) : 0.1, (5,1) : 0.2}
        params = {}
        params[JointDegreeNames.JOINT_DEGREE_TYPE] = JointDegreeType.MANUAL
        params[JointDegreeNames.JDD] = jdd
        params[JointDegreeNames.MOTIF_SIZES] = [2,3]

        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)
        
        n_vertices : int = NETWORK_SIZE
        jds = DegreeDist.sample_jds_from_jdd(n_vertices)

        params = {}
        params[GCMAlgorithmNames.MOTIF_SIZES] = [2,3]
        params[GCMAlgorithmNames.EDGE_NAMES] = ['2-clique','3-clique']
        params[GCMAlgorithmNames.BUILD_FUNCTIONS] = [clique_motif,clique_motif]
        g = GCMAlgorithmNetwork(
                params
        ).random_clustered_graph(jds)

        jdd = JointDegreeDistributionFromNetwork.get_joint_degree_distribution(g._G)

        for key in jdd.keys():
            self.assertAlmostEqual(jdd[key],jdd[key],2)
