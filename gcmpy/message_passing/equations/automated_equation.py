import itertools
import networkx as nx


class AutomatedEquation():

    def __init__(self):
        self._edge_combinations = {}
        self._connected_subgraphs = {}

    def get_connected_subgraphs(self, G: nx.Graph, root: int):
        """
        Cache for connected subgraphs of G containing root.
        """
        key: str = f"{root}-{G.name}"
        if key in self._connected_subgraphs:
            return self._connected_subgraphs[key]
        else:
            results = []
            self._get_connected_subgraphs(
                G, {root}, set(G.neighbors(root)), {root}, results, len(G.nodes())
            )
            self._connected_subgraphs[key] = results
            return self._connected_subgraphs[key]

    def _get_connected_subgraphs(self, G: nx.Graph, subgraph: set, possible: set, excluded: set, results: list, max_size: int):
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
            self._get_connected_subgraphs(G, new_subgraph, new_possible, excluded, results, max_size)


    def get_edge_combinations(self, G: nx.Graph, c: list[int]) -> list[tuple]:
        """
        Returns the edge combinations of a subgraph can be removed and still
        yield a connected subgraph. This routine is expensive, so we cache
        the combinations under the motif ID.
        """
        key: str = f"{c}-{G.name}"
        if key in self._edge_combinations:
            return self._edge_combinations[key]
        else:
            edge_combinations_final = []
            all_edge_combinations = []
            for l in range(0, len(G.edges())+1):
                all_edge_combinations.extend([e for e in itertools.combinations(G.edges(), l)])

            for es in all_edge_combinations:
                g_test: nx.Graph = G.copy()
                g_test.remove_edges_from(es)
                if nx.is_connected(g_test):
                    edge_combinations_final.append(es)

            self._edge_combinations[key] = edge_combinations_final
            return self._edge_combinations[key]


    def automated_equation(self, G: nx.Graph, p: float, root: int) -> float:
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
        components = self.get_connected_subgraphs(G, root)
        
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
            edge_combinations = self.get_edge_combinations(g, c)
            for es in edge_combinations:
                g_test = g.copy()
                g_test.remove_edges_from(es)
                us = get_us(g_test, root)
                prob += (
                    (pow(p, len(g_test.edges())) * pow(1 - p, len(es)))
                    * interface_edges
                    * us
                )

        return prob
