from gcmpy.network.network import Network
from gcmpy.gcm_algorithm.gcm_algorithm import GCMAlgorithm
from gcmpy.gcm_algorithm.gcm_algorithm_fast import FastGCMAlgorithm
from gcmpy.network.edge_list_to_network import ConvertEdgeListToNetwork
from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames


class GCMAlgorithmNetwork(GCMAlgorithm):
    def random_clustered_graph(self, jds: list) -> Network:
        params = {}
        params[GCMAlgorithmNames.MOTIF_SIZES] = self._motif_sizes
        params[GCMAlgorithmNames.BUILD_FUNCTIONS] = self._build_functions
        params[GCMAlgorithmNames.EDGE_NAMES] = self._edge_names
        CEdgeList = FastGCMAlgorithm(params).random_clustered_graph(jds)

        return ConvertEdgeListToNetwork.convert(CEdgeList)
