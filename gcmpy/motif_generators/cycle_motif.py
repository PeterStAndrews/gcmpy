from itertools import tee


def cycle_motif(vertices: list) -> list:
    """
    Accepts a list of ints and creates a cycle of adjacent vertex pairs. Creates two
    iterators and then increments one before zipping together to create tuples. Finally,
    connects the start and end of the chain toegher.
    :param vertices: list of vertices (int)
    :returns: edge list as list of tuples (int, int)
    """

    a, b = tee(vertices)
    _ = next(b, None)
    edges = list(zip(a, b))
    edges.append((vertices[0], vertices[-1]))
    return edges
