from gcmpy.network.edge_list import LightWeightEdgeList
from gcmpy.network.network import Network
from gcmpy.names.network_names import NetworkNames


class NetworkToEdgeList:
    @staticmethod
    def convert(network: Network) -> LightWeightEdgeList:
        model = LightWeightEdgeList()
        model.edge_list = list(network.G.edges())
        model.joint_degrees = [
            network.G.nodes[n][NetworkNames.JOINT_DEGREE]
            for n in range(len(network.G.nodes()))
        ]
        model.topologies = [
            network.G.edges[e][NetworkNames.TOPOLOGY] for e in network.G.edges()
        ]
        model.motif_id = [
            network.G.edges[e][NetworkNames.MOTIF_IDS] for e in network.G.edges()
        ]
        return model
