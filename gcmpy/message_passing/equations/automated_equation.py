import itertools
import networkx as nx


def get_connected_subgraphs(G: nx.Graph, subgraph: set, possible: set, excluded: set, results: list, max_size: int):
    """
    Backtracking algorithm to find all subgraphs that contain a root vertex
    from a graph represented as an edge set.
    """
    results.append(subgraph)

    if len(subgraph) == max_size:
        return

    for j in possible - excluded:
        new_subgraph: set = subgraph | {j}
        excluded: set = excluded | {j}
        new_possible: set = (possible | set(G.neighbors(j))) - excluded
        get_connected_subgraphs(G, new_subgraph, new_possible, excluded, results, max_size)


def automated_equation(G: nx.Graph, p: float, root: int) -> float:
    """
    Returns the probability that a vertex does not connect to the giant component via motif.

    Args:
        `G` (nx.Graph): is the motif.
        `p` (float): is the bond occupation prob.
        `root` (int): is the focal vertex.

    Returns:
        `prob` (float): probability `root` not in giant component via `G`
    """

    def get_us(G: nx.Graph, root: int) -> float:
        """
        Returns product of u values for each vertex in ns. Assumes that
        G has a attribute keyed by 'u', whose expected value is a float.
        """
        product = 1.0
        for n in G.nodes():
            if n == root:
                continue
            product *= G.nodes[n]["u"]
        return product

    prob = 0.0

    # `components` is a list of combinations of vertices that contain `root`
    # and that are valid subgraphs of a motif.
    components = []
    get_connected_subgraphs(G, {root}, set(G.neighbors(root)), {root}, components, len(G.nodes()))
    
    for c in components:
        
        c = list(c)
        
        if len(c) == 1:
            # isolated vertex, set all of its edges to (1-p)
            prob += pow(1 - p, len(list(G.neighbors(c[0]))))
            continue

        interface_edges = 1.0
        g = G.copy()

        # remove interface and inconsequential edges
        for e in G.edges():
            if (e[0] not in c) and (e[1] not in c):
                # remove edges between vertex pairs that are not in the same component `root`,
                # it does not contribute to the equation
                g.remove_edge(*e)
            elif (e[0] in c) and (e[1] in c):
                # edge been vertices in same component as `root`
                continue
            else:
                # one end is in component, other is out - interface edge, must be (1-p)
                g.remove_edge(*e)
                interface_edges *= (1 - p)
                continue

        # remove isolated vertices from RG
        g.remove_nodes_from([n for n in g.nodes() if len(list(g.neighbors(n))) == 0])

        # remaining edges in g are between vertices in the same component as `root`
        # they can be either p or 1-p state as long as the component is connected
        edge_combinations = []
        for l in range(0, len(g.edges())+1):
            edge_combinations.extend([e for e in itertools.combinations(g.edges(), l)])

        for es in edge_combinations:
            g_test = g.copy()
            g_test.remove_edges_from(es)

            if nx.is_connected(g_test):
                # valid edge configuration, add term to prob
                us = get_us(g_test, root)
                prob += (
                    (pow(p, len(g_test.edges())) * pow(1 - p, len(es)))
                    * interface_edges
                    * us
                )

    return prob
