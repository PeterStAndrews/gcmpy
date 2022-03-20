
from gcmpy.network.edge_list import LightWeightEdgeList
from gcmpy.network.network import Network

class NetworkToEdgeList:

    @staticmethod
    def convert(network: Network) -> LightWeightEdgeList:
        model = LightWeightEdgeList()
        model.edge_list = list(network.G.edges(),data=False)
        model.joint_degrees = [network.G.nodes[n]['joint_degree'] for n in network.G.nodes()]
        model.topologies = [network.G.edges[e]['topology'] for e in network.G.edges()]
        return model