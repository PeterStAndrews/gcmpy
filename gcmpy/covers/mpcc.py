import networkx as nx
import itertools


def MPCC(G: nx.Graph, max_size: int = 0):
    '''
    An edge-disjoint clique cover that preserves larger motifs within
    a network. This works by assigning labels to each edge of the form

    G.edges[e[0],e[1]]['clique'] = f'{len(c)}-{c}-{ID}'

    where c is a list of vertices in the clique and ID is a unique ID
    for the motif. For large networks the library call to
    `nx.enumerate_all_cliques` might be slow; however, we only need to
    call it once.

    Take care of self-loops in the network, as these will not be labelled.

    `max_size` is an optional integer to disregard cliques larger than this size.
    If it is set to zero (default) then all clique sizes are included. If it is
    greater than zero, only cliques larger than this value are included in the
    cover.

    :param G nx.Graph: the network to cover.
    :param max_size int: optional int to control the maximum clique size allowed

    :returns: G a covered graph.
    '''
    g: nx.Graph = G.copy()
    cliques: list = sorted(list(nx.enumerate_all_cliques(g)), key=len, reverse=True)
    cover: list = []
    for c in cliques:

        if len(c) > max_size and max_size > 0:
            continue

        skip: bool = False
        for e in itertools.combinations(c, 2):
            if not g.has_edge(e[0], e[1]):
                skip = True
                break
        if not skip:
            cover.append(c)
            g.remove_edges_from(list(itertools.combinations(c, 2)))

    clique_ID: int = itertools.count(0)
    for c in cover:
        ID: int = next(clique_ID)
        for e in itertools.combinations(c, 2):
            G.edges[e[0], e[1]]['clique'] = f'{len(c)}-{c}-{ID}'

    return G
