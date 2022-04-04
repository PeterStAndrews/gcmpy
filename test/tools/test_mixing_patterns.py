import unittest
from gcmpy.tools.joint_excess_joint_degree import JointExcessJointDegree
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import (
    JointDegreeManual,
)
from gcmpy.motif_generators.clique_motif import clique_motif
from gcmpy.gcm_algorithm.gcm_algorithm_network import GCMAlgorithmNetwork
from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames
from gcmpy.names.joint_degree_names import JointDegreeNames
from gcmpy.names.tools_names import ToolsNames
from gcmpy.tools.joint_excess_joint_degree_matrices import (
    JointExcessJointDegreeMatrices,
)

NETWORK_SIZE: int = 100000


class MotifMixingPatternsTest(unittest.TestCase):
    def test_two_topologies(self):

        """
        Create a network from a given degree distribution.
        Extract its ejk mixing patterns.
        Compare those to theoretical.
        """

        jdd = {(5, 1): 1 / 3, (3, 2): 1 / 3, (1, 3): 1 / 3}
        params = {}
        params[JointDegreeNames.JDD] = jdd
        params[JointDegreeNames.MOTIF_SIZES] = [2, 3]

        DegreeDistObj = JointDegreeManual(params)
        n_vertices: int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params[GCMAlgorithmNames.MOTIF_SIZES] = [2, 3]
        params[GCMAlgorithmNames.EDGE_NAMES] = ["2-clique", "3-clique"]
        params[GCMAlgorithmNames.BUILD_FUNCTIONS] = [clique_motif, clique_motif]
        g = GCMAlgorithmNetwork(params).random_clustered_graph(jds)

        params = {}
        params[ToolsNames.NETWORK] = g._G
        params[ToolsNames.EDGE_NAMES] = ["2-clique", "3-clique"]
        C = JointExcessJointDegree(params)

        ejks: JointExcessJointDegreeMatrices = C.get_ejks()

        # theoretical ejk matrices for tree and triangle
        ejk_tree = {
            (0, 3, 0, 3): 1 / 81,
            (0, 3, 4, 1): 5 / 81,
            (0, 3, 2, 2): 3 / 81,
            (4, 1, 0, 3): 5 / 81,
            (4, 1, 4, 1): 25 / 81,
            (4, 1, 2, 2): 15 / 81,
            (2, 2, 0, 3): 3 / 81,
            (2, 2, 4, 1): 15 / 81,
            (2, 2, 2, 2): 9 / 81,
        }

        ejk_triangle = {
            (3, 1, 3, 1): 16 / 144,
            (3, 1, 1, 2): 24 / 144,
            (3, 1, 5, 0): 8 / 144,
            (1, 2, 3, 1): 24 / 144,
            (1, 2, 1, 2): 36 / 144,
            (1, 2, 5, 0): 12 / 144,
            (5, 0, 3, 1): 8 / 144,
            (5, 0, 1, 2): 12 / 144,
            (5, 0, 5, 0): 4 / 144,
        }

        # check the experimental dict contains the expected keys
        # it may contain additional keys due to handshaking lemma
        # in gcm algorithm.
        for key in ejk_tree.keys():
            self.assertTrue(key in ejks.ejks["2-clique"])

        for key in ejk_triangle.keys():
            self.assertTrue(key in ejks.ejks["3-clique"])

        # check the values are similar
        for key in ejk_tree.keys():
            if key in ejks.ejks["2-clique"]:
                self.assertAlmostEqual(ejks.ejks["2-clique"][key], ejk_tree[key], 2)

        for key in ejk_triangle.keys():
            if key in ejks.ejks["3-clique"]:
                self.assertAlmostEqual(ejks.ejks["3-clique"][key], ejk_triangle[key], 2)

    def test_three_topologies(self):

        params = {}
        params[JointDegreeNames.JDD] = {
            (1, 0, 0): 0.2,
            (1, 1, 1): 0.5,
            (3, 0, 1): 0.1,
            (2, 1, 0): 0.2,
        }
        params[JointDegreeNames.MOTIF_SIZES] = [2, 3, 2]

        DegreeDistObj = JointDegreeManual(params)
        n_vertices: int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params[GCMAlgorithmNames.EDGE_NAMES] = [
            "2-clique-blue",
            "3-clique",
            "2-clique-red",
        ]
        params[GCMAlgorithmNames.BUILD_FUNCTIONS] = [
            clique_motif,
            clique_motif,
            clique_motif,
        ]
        params[GCMAlgorithmNames.MOTIF_SIZES] = [2, 3, 2]
        g = GCMAlgorithmNetwork(params).random_clustered_graph(jds)

        params = {}
        params[ToolsNames.NETWORK] = g._G
        params[ToolsNames.EDGE_NAMES] = ["2-clique-blue", "3-clique", "2-clique-red"]
        C = JointExcessJointDegree(params)

        ejk: JointExcessJointDegreeMatrices = C.get_ejks()
        self.assertTrue(len(ejk.ejks) == 3)
