import unittest
import networkx as nx

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
from gcmpy.tools.markov_chain_monte_carlo_rewiring import MarkovChainMonteCarloRewiring
from gcmpy.tools.joint_excess_from_ejk import JointExcessFromEjk
from gcmpy.tools.joint_degree_from_excess import JointDegreeFromExcess
from gcmpy.tools.joint_excess_joint_degree import JointExcessJointDegree


NETWORK_SIZE: int = 4000


class MCMCTest(unittest.TestCase):
    def test_MCMC(self):

        # for a given ejk and corresponding degree distribution,
        # create a network, rewire it and then extract its mixing
        # matrices until it is converged
        # create assorted matrices that have an equivalence class on the excess
        # joint degree and the joint degree distribution to the neutral case

        # global parameters
        edge_names: list = ["2-clique", "3-clique"]
        motif_sizes: list = [2, 3]

        epsilon_1: float = 1e-8
        epsilon_2: float = 1e-8
        epsilon_3: float = 1e-8

        ejk_tree_assorted: dict = {
            (0, 3, 0, 3): 9 / 81 - epsilon_1 - epsilon_2,
            (0, 3, 4, 1): epsilon_1,
            (0, 3, 2, 2): epsilon_2,
            (4, 1, 0, 3): epsilon_1,
            (4, 1, 4, 1): 45 / 81 - epsilon_1 - epsilon_3,
            (4, 1, 2, 2): epsilon_3,
            (2, 2, 0, 3): epsilon_2,
            (2, 2, 4, 1): epsilon_3,
            (2, 2, 2, 2): 27 / 81 - epsilon_2 - epsilon_3,
        }

        phi_1: float = 1e-8
        phi_2: float = 1e-8
        phi_3: float = 1e-8

        ejk_triangle_assorted: dict = {
            (3, 1, 3, 1): 48 / 144 - phi_1 - phi_2,
            (3, 1, 1, 2): phi_1,
            (3, 1, 5, 0): phi_2,
            (1, 2, 3, 1): phi_1,
            (1, 2, 1, 2): 72 / 144 - phi_1 - phi_3,
            (1, 2, 5, 0): phi_3,
            (5, 0, 3, 1): phi_2,
            (5, 0, 1, 2): phi_3,
            (5, 0, 5, 0): 24 / 144 - phi_2 - phi_3,
        }

        params: dict = {}
        params[ToolsNames.EDGE_NAMES] = edge_names
        params[ToolsNames.EJKS] = {
            "2-clique": ejk_tree_assorted,
            "3-clique": ejk_triangle_assorted,
        }
        ejk_target = JointExcessJointDegreeMatrices(params)

        # get excess joint degree distributions
        qks: dict = JointExcessFromEjk.get_excess_joint_distributions(ejk_target)

        # get joint degree distribution
        jdd: dict = JointDegreeFromExcess.get_joint_degree_distribution(qks, edge_names)

        # use joint degree distribution to create network
        params: dict = {}
        params[JointDegreeNames.JDD] = jdd
        params[JointDegreeNames.MOTIF_SIZES] = motif_sizes
        DegreeDistObj = JointDegreeManual(params)
        jds: list = DegreeDistObj.sample_jds_from_jdd(NETWORK_SIZE)

        params: dict = {}
        params[GCMAlgorithmNames.MOTIF_SIZES] = motif_sizes
        params[GCMAlgorithmNames.EDGE_NAMES] = edge_names
        params[GCMAlgorithmNames.BUILD_FUNCTIONS] = [clique_motif, clique_motif]
        g = GCMAlgorithmNetwork(params).random_clustered_graph(jds)

        # pull out the randomly mixed matrices
        params: dict = {}
        params[ToolsNames.NETWORK] = g.G
        params[ToolsNames.EDGE_NAMES] = edge_names
        C = JointExcessJointDegree(params)
        initial_experimental_ejks: JointExcessJointDegreeMatrices = C.get_ejks()

        # rewire the network to target via MCMC
        params: dict = {}
        params[ToolsNames.NETWORK] = g
        params[ToolsNames.EJKS] = ejk_target
        params[ToolsNames.SEARCH_LIMIT] = 20
        params[ToolsNames.CONVERGENCE_LIMIT] = 25000
        mcmc = MarkovChainMonteCarloRewiring(params)
        G: nx.Graph = mcmc.rewire()

        # extract the joint excess joint degree distribution from the graph
        params: dict = {}
        params[ToolsNames.NETWORK] = G
        params[ToolsNames.EDGE_NAMES] = edge_names
        CJointExcess = JointExcessJointDegree(params)
        experimental_ejks: JointExcessJointDegreeMatrices = CJointExcess.get_ejks()

        # check that the experimental mixing matrices are closer to the target than the initial
        # check that the qks and jdd are correct to theoretical

        theoretical_jdd: dict = {
            (5, 1): 1 / 3, (3, 2): 1 / 3, (1, 3): 1 / 3
        }

        theoretical_qk_tree: dict = {
            (0, 3): 0.11025028370284273,
            (4, 1): 0.5568150524081263,
            (2, 2): 0.33293067275708305,
        }

        theoretical_qk_triangle: dict = {
            (3, 1): 0.33423484065484,
            (1, 2): 0.4980622601010759,
            (5, 0): 0.16769688919908976,
        }

        self.assertTrue(jdd == theoretical_jdd)

        for key in qks["2-clique"]:
            self.assertAlmostEqual(qks["2-clique"][key], theoretical_qk_tree[key], 2)

        for key in qks["3-clique"]:
            self.assertAlmostEqual(
                qks["3-clique"][key], theoretical_qk_triangle[key], 2
            )

        # check mixing matrices
        print('target', 'initial', 'final')
        for topology in edge_names:
            for key in ejk_target.ejks[topology]:
                target: float = ejk_target.ejks[topology][key]
                initial: float = initial_experimental_ejks.ejks[topology][key]
                final: float = experimental_ejks.ejks[topology][key]
                self.assertTrue(abs(target - initial) >= abs(target - final))
                print(target, initial, final)
