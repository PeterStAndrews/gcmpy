import unittest
from gcmpy.gcm_algorithm.gcm_algorithm_custom_motifs import GCMAlgorithmCustomMotifs
from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames


class GCMAlgorithmCustomMotifsTest(unittest.TestCase):
    def test_custom_motifs(self):
        def diamond(vs):
            return (
                (vs[0], vs[1]),
                (vs[1], vs[2]),
                (vs[2], vs[3]),
                (vs[3], vs[1]),
                (vs[0], vs[2]),
            )

        def diamond_names():
            return (
                "diamond-outer",
                "diamond-outer",
                "diamond-outer",
                "diamond-outer",
                "diamond-inner",
            )

        def twoclique(vs):
            return (vs[0], vs[1])

        def twoclique_names():
            return "2-clique"

        def threeclique(vs):
            return (vs[0], vs[1]), (vs[0], vs[2]), (vs[1], vs[2])

        def threeclique_names():
            return "3-clique", "3-clique", "3-clique"

        def pentagon_in_manuscript(vs):
            return (
                (vs[0], vs[1]),
                (vs[1], vs[2]),
                (vs[2], vs[3]),
                (vs[3], vs[4]),
                (vs[0], vs[4]),
                (vs[1], vs[3]),
            )

        def pentagon_in_manuscript_names():
            return "p01", "p12", "p23", "p34", "p40", "p13"

        # two cliques, three cliques, diamonds and pentagon in paper
        # must satisfy handshaking lemma and Erdos-Gallai inequality suitably generalised.
        jds = [
            (2, 1, 0, 1, 1, 0, 0),
            (1, 1, 0, 1, 1, 0, 0),
            (3, 1, 1, 0, 0, 1, 0),
            (2, 0, 1, 0, 0, 1, 0),
            (0, 0, 0, 1, 0, 0, 1),
            (1, 0, 0, 1, 0, 0, 0),
            (1, 0, 1, 0, 0, 0, 0),
            (1, 0, 1, 0, 0, 0, 0),
            (1, 0, 0, 1, 0, 0, 0),
            (1, 0, 0, 1, 0, 0, 0),
            (1, 0, 1, 0, 0, 0, 0),
            (0, 0, 1, 0, 0, 0, 0),
        ]

        params = {}

        # number of vertices of this orbit in motif
        params[GCMAlgorithmNames.MOTIF_SIZES] = [2, 3, 2, 2, 2, 2, 1]

        # function callbacks to return names of edges
        params[GCMAlgorithmNames.EDGE_NAMES] = [
            twoclique_names,
            threeclique_names,
            diamond_names,
            pentagon_in_manuscript_names,
        ]

        # function callbacks to return the edges
        params[GCMAlgorithmNames.BUILD_FUNCTIONS] = [
            twoclique,
            threeclique,
            diamond,
            pentagon_in_manuscript,
        ]

        # indices of joint degree tuple that belong to each motif
        # in this case, 2-clique, 3-clique, diamond and pentagon cycle with chord (1,3)
        params[GCMAlgorithmNames.MOTIF_INDICES] = [[0], [1], [2, 3], [4, 5, 6]]

        es = GCMAlgorithmCustomMotifs(params).random_clustered_graph(jds)

        self.assertTrue(len(es.edge_list) == len(es.topologies) == 31)
