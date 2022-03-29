
import unittest

from gcmpy.network.edge_list_to_network import ConvertEdgeListToNetwork
from gcmpy.network.network_to_edge_list import NetworkToEdgeList
from gcmpy.names.joint_degree_names import JointDegreeNames
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import JointDegreeManual
from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames
from gcmpy.gcm_algorithm.gcm_algorithm_network import GCMAlgorithmNetwork
from gcmpy.motif_generators.clique_motif import clique_motif
from gcmpy.network.edge_list import LightWeightEdgeList
from gcmpy.network.network import Network
from gcmpy.names.network_names import NetworkNames


NETWORK_SIZE = 100000

class ConversionTests(unittest.TestCase):

    def test_conversions(self):

        # create a random network
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

        # convert it to an edge list
        edge_list: LightWeightEdgeList = NetworkToEdgeList.convert(g)

        # convert it back again to a network
        g2: Network = ConvertEdgeListToNetwork.convert(edge_list)

        # test equality
        for e in g._G.edges():
            self.assertTrue(g2.G.has_edge(*e))
            self.assertTrue(
                g.G.edges[e][NetworkNames.TOPOLOGY] == g2.G.edges[e][NetworkNames.TOPOLOGY]
            )
            self.assertTrue(
                g.G.edges[e][NetworkNames.MOTIF_IDS] == g2.G.edges[e][NetworkNames.MOTIF_IDS]
            )

        for n in g.G.nodes():
            self.assertTrue(
                g.G.nodes[n][NetworkNames.JOINT_DEGREE] == g2.G.nodes[n][NetworkNames.JOINT_DEGREE]
            )
        
        