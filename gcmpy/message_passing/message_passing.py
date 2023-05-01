import networkx as nx
from gcmpy.message_passing.message_passing_mixin import MessagePassingMixin
from gcmpy.message_passing.equations.automated_equation import AutomatedEquation


class MessagePassing:
    def __init__(
        self, G: nx.Graph, cover_type: str = 'motif cover', iterations: int = 25
    ):
        """
        An implementation of the message passing algorithm for networkx graphs that
        have been covered with edge-disjoint motifs.
        
        :param G: networkx graph with edge labels.
        :param cover_type: [optional] str, the type of cover being calculated
        :param iterations: [optional] integer number of iterations for fixed point calculation
        """
        self._MPM = MessagePassingMixin(cover_type, G)
        self._AE = AutomatedEquation()
        self._H_tau: dict = {}
        self._iterations: int = iterations

    def resolve_equation(self, focal: int, label: str, prods: dict) -> float:
        """
        Calculates the probability that connection to the GCC fails through this motif. The
        implementation of this method depends on what is being calculated. If the cover only
        contains cliques, then we can do the following

        ```python
            topology = self._MPM.get_motif_topology(label)
            return clique_equation(topology, self._phi, prods.values())
        ```
    
        Another plan would be to store a list of callbacks, which could be provided on
        initialisation of this class

        ```python
            return self._equations[topology](args)
        ```

        A final option is to call the automated equation, which will require a nx subgraph
        to be constructed first and the u = prod_{\nu}H_{j\leftarrow \nu} values installed 
        on the vertices as attributes keyed by `u`.
        
        ```python
            edges: list[tuple] = self._MPM.get_edges_in_motif(label)
            H = nx.Graph()
            H.add_edges_from(edges)
            nx.set_node_attributes(H, prods, "u")
            return automated_equation(H, self._phi, focal)
        ```
    
        :param focal: int, the ID of the focal vertex
        :param label: str, cover label for an edge in the motif
        :param prods: a dict of floats of `j: H_{j leftarrow nu}(z)'
        :return float: the probability that connection to the GCC fails through this motif.
        """
        edges: list[tuple] = self._MPM.get_edges_in_motif(label)
        H = nx.Graph(name=f"{focal}-{self._MPM.get_motif_ID(label)}")
        H.add_edges_from(edges)
        nx.set_node_attributes(H, prods, "u")
        return self._AE.automated_equation(H, self._phi, focal)
        
    def calculate_H_tau(
        self, focal: int, label: str
    ) -> None:
        """
        Calculates `H_{focal leftarrow tau}(z)' which is the probability that `focal' vertex
        does not become attached to the GCC from membership in motif tau. 
        `vertices_in_motif` are vertices that belong to *this* motif.

        :param focal: int the ID of the focal vertex
        :param label: str of the cover label for this motif
        """
        # parse the topology of the motif, its unique ID and
        # the other vertices in the motif.
        motif_ID: str = self._MPM.get_motif_ID(label)
        vertices_in_motif: tuple = self._MPM.get_vertices_in_motif(label)

        # iterate the vertices of *this* motif apart from focal vertex i
        prods = {}
        for j in vertices_in_motif:
            if j == focal:
                continue

            # Get all of the direct neighbours of j. Note non-direct neighbours
            # in motif nu will still be caught by cover label.
            js_neighbours = list(self._MPM._G.neighbors(j))

            # ignore verticies in *this* motif with focal
            js_neighbours = set(js_neighbours) - set(vertices_in_motif)

            # for each motif of `j' get prob don't connect j to GCC
            prod_j = 1
            done_motifs = set()
            for l in js_neighbours:
                label_l = self._MPM.get_edge_cover_label(j, l)
                motif_ID_l = self._MPM.get_motif_ID(label_l)

                if motif_ID_l in done_motifs:
                    continue

                prod_j *= self._H_tau[(j, motif_ID_l)]
                done_motifs.add(motif_ID_l)

            # H_{\tau_j\leftarrow l}(z)
            prods[j] = prod_j

        # H_{i\leftarrow \tau}(z)
        self._H_tau[(focal, motif_ID)] = self.resolve_equation(focal, label, prods)

    def theoretical(self, phi: float) -> float:
        """
        Performs the message passing algorithm for networks with an
        edge-disjoint cover of motifs for bond percolation dynamics.

        :param phi: float, bond occupation probability
        :return float: average probability that a random vertex belongs to the GCC
        """
        self._phi = phi

        # initialise the model
        self._H_tau: dict = {}
        for i, j in self._MPM._G.edges():
            label: str = self._MPM.get_edge_cover_label(i, j)
            motif_ID: str = self._MPM.get_motif_ID(label)

            for k in self._MPM.get_vertices_in_motif(label):
                self._H_tau[(k, motif_ID)] = 0.5

        # fixed point iteration repeats
        for _ in range(self._iterations):
            # iterate all edges in the network
            for i, j in self._MPM._G.edges():
                # pull the cover label from the raw network label
                label: str = self._MPM.get_edge_cover_label(i, j)

                # ================= Calculate H_{i,\tau}(z) ================== #
                self.calculate_H_tau(i, label)

                # ================= Calculate H_{j,\tau}(z) ================== #
                self.calculate_H_tau(j, label)

        # construct the probability that the average vertex does not
        # belong to the GCC.
        outer_sum = 0.0
        for i in self._MPM._G.nodes():
            prod = 1
            done_motifs = set()
            for j in self._MPM._G.neighbors(i):
                label = self._MPM.get_edge_cover_label(i, j)
                motif_ID = self._MPM.get_motif_ID(label)

                if motif_ID in done_motifs:
                    continue

                prod *= self._H_tau[(i, motif_ID)]
                done_motifs.add(motif_ID)

            outer_sum += prod

        return 1 - ((1.0 * outer_sum) / self._MPM._G.order())
