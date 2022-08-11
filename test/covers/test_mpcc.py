import unittest
import networkx as nx

from gcmpy.covers.mpcc import MPCC
from gcmpy.motif_generators.clique_motif import clique_motif
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import (
    JointDegreeManual,
)
from gcmpy.names.joint_degree_names import JointDegreeNames
from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames
from gcmpy.gcm_algorithm.gcm_algorithm_network import GCMAlgorithmNetwork


NETWORK_SIZE: int = 100000


class MPCCTest(unittest.TestCase):

    def test_mpcc(self):

        params = {}
        params[JointDegreeNames.JDD] = {
            (1, 0): 0.2,
            (2, 1): 0.5,
            (3, 0): 0.1,
            (5, 1): 0.2,
        }
        params[JointDegreeNames.MOTIF_SIZES] = [2, 3]

        DegreeDistObj = JointDegreeManual(params)
        n_vertices: int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        params = {}
        params[GCMAlgorithmNames.MOTIF_SIZES] = [2, 3]
        params[GCMAlgorithmNames.EDGE_NAMES] = ["2-clique", "3-clique"]
        params[GCMAlgorithmNames.BUILD_FUNCTIONS] = [clique_motif, clique_motif]
        g = GCMAlgorithmNetwork(params).random_clustered_graph(jds)

        # the MPCC method does not work with self loops, so this is required!
        g.G.remove_edges_from(nx.selfloop_edges(g.G))

        # Extract the MPCC clique cover (all size cliques allowed)
        g.G = MPCC(g.G)

        # iterate the edges and record the clique size
        cliques = {}
        for e in g.G.edges():
            label = g.G.edges[e[0], e[1]]['clique']
            size = int(label.split('-')[0])
            cliques[size] = cliques.get(size, 0) + 1

        for size in cliques:
            cliques[size] /= size

        self.assertTrue(max([size for size in cliques]) == 3)

        # Extract the 2-clique cover
        G = MPCC(g.G, 2)

        # iterate the edges and record the clique size
        cliques = {}
        for e in G.edges():
            label = G.edges[e[0], e[1]]['clique']
            size = int(label.split('-')[0])
            cliques[size] = cliques.get(size, 0) + 1

        for size in cliques:
            cliques[size] /= size

        self.assertTrue(max([size for size in cliques]) == 2)
