
from itertools import combinations

def clique_motif(nodes: list) -> list:
    """
    Accepts list of ints and creates all possible pairs which it returns as a list of tuples of int pairs
    :param nodes: list of nodes (int)
    :returns: edge list as list of tuples (int, int)
    """
    return list(combinations(nodes,2))