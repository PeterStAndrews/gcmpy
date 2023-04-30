import unittest
import networkx as nx

from gcmpy.message_passing.number_connected_graphs import number_of_connected_graphs


class NumberConnectedGraphsTest(unittest.TestCase):

    def test_number_connected_graphs(self):

        G = nx.Graph()

        G.add_edge(0, 1)
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(3, 4)
        G.add_edge(4, 0)

        G.add_edge(2, 0)
        G.add_edge(2, 4)
        G.add_edge(3, 1)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 1
        self.assertTrue(3 == number_of_connected_graphs(substrate_motif,
                                                        neighbours_in_component,
                                                        focal_vertex,
                                                        number_of_edges_to_remove)
                        )

        G = nx.complete_graph(6)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2, 3, 4, 5]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 10

        num: int = number_of_connected_graphs(substrate_motif,
                                              neighbours_in_component,
                                              focal_vertex,
                                              number_of_edges_to_remove)

        self.assertTrue(1296 == num)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2, 3, 4, 5]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 9

        num: int = number_of_connected_graphs(substrate_motif,
                                              neighbours_in_component,
                                              focal_vertex,
                                              number_of_edges_to_remove)

        self.assertTrue(3660 == num)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2, 3, 4, 5]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 1

        num: int = number_of_connected_graphs(substrate_motif,
                                              neighbours_in_component,
                                              focal_vertex,
                                              number_of_edges_to_remove)

        self.assertTrue(15 == num)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2, 3, 4]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 4

        num: int = number_of_connected_graphs(substrate_motif,
                                              neighbours_in_component,
                                              focal_vertex,
                                              number_of_edges_to_remove)

        self.assertTrue(205 == num)


class AutomatedNumberConnectedGraphs(unittest.TestCase):

    def test_number_connected_graphs(self):

        G = nx.Graph()

        G.add_edge(0, 1)
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(3, 4)
        G.add_edge(4, 0)

        G.add_edge(2, 0)
        G.add_edge(2, 4)
        G.add_edge(3, 1)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 1
        self.assertTrue(3 == number_of_connected_graphs(substrate_motif,
                                                        neighbours_in_component,
                                                        focal_vertex,
                                                        number_of_edges_to_remove)
                        )

        G = nx.complete_graph(6)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2, 3, 4, 5]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 10

        num: int = number_of_connected_graphs(substrate_motif,
                                              neighbours_in_component,
                                              focal_vertex,
                                              number_of_edges_to_remove)

        self.assertTrue(1296 == num)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2, 3, 4, 5]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 9

        num: int = number_of_connected_graphs(substrate_motif,
                                              neighbours_in_component,
                                              focal_vertex,
                                              number_of_edges_to_remove)

        self.assertTrue(3660 == num)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2, 3, 4, 5]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 1

        num: int = number_of_connected_graphs(substrate_motif,
                                              neighbours_in_component,
                                              focal_vertex,
                                              number_of_edges_to_remove)

        self.assertTrue(15 == num)

        substrate_motif: nx.Graph = G
        neighbours_in_component: list = [1, 2, 3, 4]
        focal_vertex: int = 0
        number_of_edges_to_remove: int = 4

        num: int = number_of_connected_graphs(substrate_motif,
                                              neighbours_in_component,
                                              focal_vertex,
                                              number_of_edges_to_remove)

        self.assertTrue(205 == num)