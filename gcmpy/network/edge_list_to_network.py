import networkx as nx

from gcmpy.network.edge_list import LightWeightEdgeList
from gcmpy.network.network import Network
from gcmpy.names.network_names import NetworkNames


class EdgeListToNetwork:
    @staticmethod
    def convert(edgelist: LightWeightEdgeList) -> Network:

        model = Network()
        model.G.add_edges_from(edgelist.edge_list)

        # create vertex attributes dict
        joint_degrees = {}
        for n, jd in enumerate(edgelist.joint_degrees):
            joint_degrees[n] = jd
        nx.set_node_attributes(model.G, joint_degrees, NetworkNames.JOINT_DEGREE)
        # create edge attributes dict
        topologies = {}
        motif_ids = {}
        for e, name, motif_id in zip(
            edgelist.edge_list, edgelist.topologies, edgelist.motif_id
        ):
            topologies[e] = name
            motif_ids[e] = motif_id

        nx.set_edge_attributes(model.G, topologies, NetworkNames.TOPOLOGY)
        nx.set_edge_attributes(model.G, motif_ids, NetworkNames.MOTIF_IDS)

        return model
