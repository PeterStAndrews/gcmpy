


from gcmpy.network.network import Network
from gcmpy.gcm_algorithm.gcm_algorithm import GCMAlgorithm
from gcmpy.gcm_algorithm.gcm_algorithm_fast import FastGCMAlgorithm
from gcmpy.network.edge_list_to_network import ConvertEdgeListToNetwork

class GCMAlgorithmNetwork(GCMAlgorithm):

    def random_clustered_graph(self, jds: list) -> Network:
        params = {}
        params["motif_sizes"] = self._motif_sizes
        params["build_functions"] = self._build_functions
        params["edge_names"] = self._edge_names
        CEdgeList = FastGCMAlgorithm(
            params
        ).random_clustered_graph(jds)

        return ConvertEdgeListToNetwork.convert(CEdgeList)