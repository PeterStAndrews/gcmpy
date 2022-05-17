import random
import networkx as nx

def bond_percolate(g: nx.Graph, phi: float) -> float:
    """
    Performs bond percolation with occupancy probability `phi`
    :param g: nx graph to percolate (original)
    :param phi: bond occupancy probability
    :returns float S: size of the GCC
    """
    G: nx.Graph = g.copy()
    es: list = [e for e in G.edges() if random.random() > phi]
    G.remove_edges_from(es)
    Gcc: list = sorted(nx.connected_components(G), key=len, reverse=True)
    return float(len(Gcc[0]))/G.order()