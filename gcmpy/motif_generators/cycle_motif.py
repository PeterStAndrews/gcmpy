
from itertools import tee

def cycle_motif(nodes: list) -> list:
    """
    Accepts a list of ints and creates a cycle of adjacent node pairs. Creates two 
    iterators and then increments one before zipping together to create tuples. Finally,
    connects the start and end of the chain toegher.
    :param nodes: list of nodes (int)
    :returns: edge list as list of tuples (int, int)
    """
    
    a, b = tee(nodes)
    _ = next(b, None)
    edges = list(zip(a,b))
    edges.append((nodes[0], nodes[-1]))
    return edges