import networkx as nx


class Network:
    def __init__(self):
        self._G = nx.Graph()

    def add_edge(self, e: tuple):
        """
        Adds an edge to the networkx
        :param e tuple: edge
        """
        self._G.add_edge(*e)

    def add_edges_from(self, edges: list) -> None:
        """
        Adds edges from list of tuples (int,int) to the edge list.
        :param edges: list of tuples of ints.
        """
        self._G.add_edges_from(edges)

    def find_cliques(self):
        """
        Returns all maximal cliques
        """
        return list(nx.find_cliques(self._G))

    def remove_edge(self, i: int, j: int) -> None:
        """
        Removes edge (i,j) from G. With nx, will also remove (j,i).
        """
        try:
            self._G.remove_edge(i, j)
        except nx.NetworkXError:
            return

    def has_edges(self) -> bool:
        """
        True if graph has edges remaining.
        """
        return len(self._G.edges()) > 0

    @property
    def G(self) -> nx.Graph:
        return self._G

    @G.setter
    def G(self, value: nx.Graph) -> None:
        self._G = value
