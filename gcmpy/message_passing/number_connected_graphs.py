import itertools
import networkx as nx


def number_of_connected_graphs(G: nx.Graph, ak: list, i: int, k: int):
    '''
    Returns the number of connected graphs among vertices
    (ak + i) with k edges removed. For instance, there are
    three graphs among a triangle with one edge removed.

    :param G: the substrate subgraph
    :param ak: the vertices in G that i connects to
    :param i: the focal vertex
    :param k: the number of edges to remove

    :returns the number of connected graphs that contain vertices
    ak+i, with k edges removed.
    '''

    count: int = 0
    H: nx.Graph = G.copy()
    for n in G.nodes():
        if n == i or n in ak:
            continue

        H.remove_node(n)

    for comb in itertools.combinations(H.edges(), k):
        J: nx.Graph = H.copy()
        for e in comb:
            J.remove_edge(*e)
        if nx.is_connected(J):
            count += 1

    return count
