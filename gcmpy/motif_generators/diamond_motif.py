from gcmpy.motif_generators.cycle_motif import cycle_motif


def diamond_motif(vertices: list) -> list:
    """
    Accepts list of ints and creates edge pairs for a diamond motif of 4 vertices.
    :param vertices: list of vertices (int)
    :returns: edge list as list of tuples (int, int)
    """
    if len(vertices) < 4:
        raise ("Error during motif construction - diamond_motif")

    n0, n1, n2, n3 = vertices
    edges = cycle_motif(vertices)
    edges.append((n0, n2))
    edges.append((n1, n3))
    return edges
