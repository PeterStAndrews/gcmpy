

import networkx as nx

from gcmpy.network.edge_list import LightWeightEdgeList
from gcmpy.network.network import Network

class ConvertEdgeListToNetwork:
    @staticmethod
    def convert(edgelist: LightWeightEdgeList) -> Network:

        model = Network()
        model.G.add_edges_from(edgelist.edge_list)

        # create node attributes dict
        joint_degrees = {}
        for n, jd in enumerate(edgelist.joint_degrees):
            joint_degrees[n] = jd
        nx.set_node_attributes(
            model.G, joint_degrees, "joint_degree"
        )
        # create edge attributes dict
        topologies = {}
        for e, name in zip(edgelist.edge_list, edgelist.topologies):
            topologies[e] = name
        nx.set_edge_attributes(
            model.G, topologies,"topology"
        )

        return model




