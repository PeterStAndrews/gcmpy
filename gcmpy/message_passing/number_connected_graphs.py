import itertools
import networkx as nx
from functools import lru_cache
from math import factorial


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


@lru_cache(maxsize=None)
def QQ(n: int, k: int) -> int:
    '''
    Returns the number of connected graphs of n vertices with k edges.
    This function is for cliques only.

    This function is slow since it calculates by brute force. For a
    different algorithm, see `Q(n,k)` below.

    :param n: integer number of vertices
    :param k: integer number of edges

    :returns int: the number of connected graphs with `n` vertices and
    `k` edges.

    test suite
    ----------

    checks = [1,15,105,455,1365,2997,4945,6165,5700,3660,1296,0]
    for i in range(0,len(checks)):
        print(QQ(6,15-i), checks[i])

    '''
    gg = nx.complete_graph(n)
    focal = 0
    aks = [x for x in range(1, n)]
    all_edges = int(0.5*n*(n-1))
    edges_to_remove = all_edges - k
    return number_of_connected_graphs(gg, aks, focal, edges_to_remove)


@lru_cache(maxsize=None)
def binomial(n, k):
    d = n - k
    if d < 0:
        return 0
    return factorial(n) // factorial(k) // factorial(d)


@lru_cache(maxsize=None)
def Q(n: int, k: int) -> int:
    '''
    Number of labeled, simply connected Graphs of order n, size k.
    This is a recursion taken from "Frank Harary and Edgar M. Palmer.
    Graphical Enumeration. Academic Press. 1973".

    :param n: integer number of vertices
    :param k: integer number of edges

    :returns int: the number of connected graphs with `n` vertices and
    `k` edges.

    '''
    s = n * (n - 1) // 2
    if k < n - 1 or k > s:
        res = 0
    elif k == n - 1:
        res = int(pow(n, (n - 2)))
    else:
        res = binomial(s, k)
        for m in range(0, n - 1):
            res1 = 0
            lb = max(0, k - (m + 1) * m // 2)
            for p in range(lb, k - m + 1):
                np = (n - 1 - m) * (n - 2 - m) // 2
                res1 += binomial(np, p) * Q(m + 1, k - p)

            res -= binomial(n - 1, m) * res1

    return res
