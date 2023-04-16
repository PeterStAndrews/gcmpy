import networkx as nx


class JointExcessDegree:
    @staticmethod
    def get_ejk(G: nx.Graph) -> dict:
        """
        Returns the joint excess degree (overall) of a
        networkx graph.

        :param G: networkx graph
        :return ejk: dict {(j,k) : val}
        """

        ejk = {}

        num_edges = len(G.edges())
        for e in G.edges():
            u, v = e
            u_degree = G.degree(u)
            v_degree = G.degree(v)

            u_excess_degree = u_degree - 1
            v_excess_degree = v_degree - 1

            key1 = (u_excess_degree, v_excess_degree)
            key2 = (v_excess_degree, u_excess_degree)

            ejk[key1] = ejk.get(key1, 0.0) + (0.5 / num_edges)
            ejk[key2] = ejk.get(key2, 0.0) + (0.5 / num_edges)

        return ejk
