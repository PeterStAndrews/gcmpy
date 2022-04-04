from itertools import combinations


def clique_motif(vertices: list) -> list:
    """
    Accepts list of ints and creates all possible pairs which it returns as a list of tuples of int pairs
    :param vertices: list of vertices (int)
    :returns: edge list as list of tuples (int, int)
    """
    return list(combinations(vertices, 2))
