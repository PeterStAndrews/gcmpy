import networkx as nx

from gcmpy.message_passing.message_passing_mixin import MessagePassingMixin


class MessagePassing():

    def __init__(self, cover_type: str, G: nx.Graph, equations: dict, iterations: int = 25):
        '''
        An implementation of the message passing algorithm for networkx graphs that
        have been covered with edge-disjoint motifs.

        The equations dict is a key:value pair of callback functions that, for a given
        motif topology (key) evaluates the probability that the motif fails to connect
        a vertex to the GCC. The signature of the callbacks is fixed, for instance for
        ordinary edges we have

            def func(phi: float, Hs: list[float])->float:
                return 1 - phi + phi*Hs[0]

        :param cover_type: str, the type of cover being calculated
        :param G: networkx graph with edge labels.
        :param equations: dict of callbacks that evaluate the motif equations.
        :param iterations: [optional] integer number of iterations for fixed point calculation
        '''
        self._MPM = MessagePassingMixin(cover_type, G)
        self._H_tau: dict = {}
        self._equations: dict[callable] = equations
        self._iterations: int = iterations

    def resolve_equation(self, topology: str, prods: list) -> float:
        '''
        Calculates the probability that connection to the GCC
        fails through this motif. This will throw if the equation
        is not found for the motif topology.

        :param topology: str, the toplogy of the motif
        :param prods: a list of floats of `H_{j leftarrow nu}(z)'

        :return float: the probability that connection to the GCC
        fails through this motif.
        '''
        return self._equations[topology](self._phi, prods)

    def calculate_H_tau(self, focal: int, vertices_in_motif: list, motif_ID: int,
                        topology: str) -> None:
        '''
        Calculates `H_{focal leftarrow tau}(z)' which is the probability that `focal' vertex
        does not become attached to the GCC from membership in motif tau.

        :param focal: int the ID of the focal vertex
        :param vertices_in_motif: vertices that belong to this motif
        :param motif_ID: integer motif unique ID
        :param topology: string indicating the topology
        '''
        # iterate the vertices of *this* motif apart from focal vertex i
        prods = []
        for vertex in vertices_in_motif:

            if vertex == focal:
                continue

            # get all of the neighbours of the vertex
            cavity_neighbours = [n for n in self._MPM._G.neighbors(vertex)]

            # remove vertices that are in the motif itself
            try:
                for v in vertices_in_motif:
                    if v != vertex:
                        cavity_neighbours.remove(v)
            except Exception:
                # if we are here, we have tried to remove a neighbour
                # that is not actually directly connected to `vertex`
                # but *is* in the motif.

                # it might also be an error in the configuration model
                # due to a mis-formed motif (self-loops etc).

                # just do nothing (no-op)
                pass

            # for each neighbour of `vertex' (which is focal's neighbour in the motif)
            prod_vertex = 1
            done_motifs = set()
            for ell in cavity_neighbours:

                label_l = self._MPM.get_edge_cover_label(vertex, ell)
                motif_ID_l = self._MPM.get_motif_ID(label_l)

                if motif_ID_l in done_motifs:
                    continue

                prod_vertex *= self._H_tau[(vertex, motif_ID_l)]
                done_motifs.add(motif_ID_l)

            # H_{\tau_j\leftarrow l}(z)
            prods.append(prod_vertex)

        # H_{i\leftarrow \tau}(z)
        self._H_tau[(focal, motif_ID)] = self.resolve_equation(topology, prods)

    def theoretical(self, phi: float) -> float:
        '''
        Performs the message passing algorithm for networks with an
        edge-disjoint cover of motifs for bond percolation dynamics.

        :param phi: float, bond occupation probability
        :return float: average probability that a random vertex belongs to the GCC
        '''
        self._phi = phi

        # initialise the model
        self._H_tau: dict = {}
        for i, j in self._MPM._G.edges():

            label: str = self._MPM.get_edge_cover_label(i, j)
            motif_ID: str = self._MPM.get_motif_ID(label)

            for k in self._MPM.get_vertices_in_motif(label):
                self._H_tau[(k, motif_ID)] = 0.5

        # fixed point iteration repeats
        for r in range(self._iterations):

            # iterate all edges in the network
            for i, j in self._MPM._G.edges():

                # pull the cover label from the raw network label
                label: str = self._MPM.get_edge_cover_label(i, j)

                # parse the topology of the motif, its unique ID and
                # the other vertices in the motif.
                topology: str = self._MPM.get_motif_topology(label)
                motif_ID: str = self._MPM.get_motif_ID(label)
                vertices_in_motif: tuple = self._MPM.get_vertices_in_motif(label)

                # ================= Calculate H_{i,\tau}(z) ================== #
                self.calculate_H_tau(i, vertices_in_motif, motif_ID, topology)

                # ================= Calculate H_{j,\tau}(z) ================== #
                self.calculate_H_tau(j, vertices_in_motif, motif_ID, topology)

        # construct the probability that the average vertex does not
        # belong to the GCC.
        outer_sum = 0.0
        for i in self._MPM._G.nodes():
            prod = 1
            done_motifs = set()
            for j in self._MPM._G.neighbors(i):

                label = self._MPM.get_edge_cover_label(i, j)
                topology = self._MPM.get_motif_topology(label)
                motif_ID = self._MPM.get_motif_ID(label)

                if motif_ID in done_motifs:
                    continue

                prod *= self._H_tau[(i, motif_ID)]
                done_motifs.add(motif_ID)

            outer_sum += prod

        return 1 - ((1.0*outer_sum)/self._MPM._G.order())
