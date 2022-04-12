import random
import networkx as nx
from logging import Logger

from gcmpy.names.tools_names import ToolsNames
from gcmpy.network.network import Network
from gcmpy.names.network_names import NetworkNames
from gcmpy.tools.joint_excess_joint_degree_matrices import (
    JointExcessJointDegreeMatrices,
)
from gcmpy.tools.proposal_edge import ProposalEdge
from gcmpy.tools.markov_chain_monte_carlo import MarkovChainMonteCarlo
from gcmpy.tools.joint_excess_joint_degree_keys_view import (
    JointExcessJointDegreeKeysView,
)
from gcmpy.tools.draw_set import DrawSet

class ErrorMarkovChainMonteCarloRewiring(Exception):
    ...


class MarkovChainMonteCarloRewiring(MarkovChainMonteCarlo):

    """
    Stochastically rewires `_network` so that the target correlations match `ejks`.
    """

    def __init__(self, params: dict):
        self._network: Network = None
        self._convergence_limit: int = None
        self._search_limit: int = 25
        self._ejks: JointExcessJointDegreeMatrices = None
        self._logger: Logger = Logger("MarkovChainMonteCarloRewiring")
        self._proposal_edges: list[ProposalEdge] = []
        self._proposal_count: int = 0
        self._proposals_accepted: int = 0
        self._acceptance_ratio: list = []
        try:
            self._network: Network = params[ToolsNames.NETWORK]
            self._ejks: JointExcessJointDegreeMatrices = params[ToolsNames.EJKS]
            if ToolsNames.CONVERGENCE_LIMIT in params:
                self._convergence_limit: int = params[ToolsNames.CONVERGENCE_LIMIT]
            else:
                # set a limit for the number of swaps before termination
                self._convergence_limit = 10 * self._network.G.edges()
            if ToolsNames.SEARCH_LIMIT in params:
                self._search_limit: int = params[ToolsNames.SEARCH_LIMIT]
        except Exception as e:
            raise ErrorMarkovChainMonteCarloRewiring(
                f"Error instantiating \
                {self.__class__.__name__}: {e}"
            )

    def get_all_edges(self, G: nx.Graph, u0: int, edge: tuple) -> list:
        """
        Returns all edges that have a given motif id as list of ordered tuples
        such that e = (u0,u1).
        :param G: nx.Graph
        :param e: edge in motif
        :returns list: all edges in motif involving focal vertex u0
        """
        es = []
        motif_id: int = G.edges[edge][NetworkNames.MOTIF_IDS]
        for e in G.edges(u0):
            if G.edges[e][NetworkNames.MOTIF_IDS] == motif_id:
                es.append((u0, self.get_other_vertex(u0, e)))
        return es

    def get_hashmap(self, G: nx.Graph, es: list) -> dict[list]:
        """
        Creates a hashmap of edges keyed by their topology.
        :param G: nx graph
        :param es: list of edges
        """
        hashmap_es: dict = {}
        for e in es:
            topology: str = G.edges[e][NetworkNames.TOPOLOGY]
            lst: list = hashmap_es.get(topology, [])
            lst.append(e)
            hashmap_es[topology] = lst
        return hashmap_es

    def is_edge_choice_suitable(
        self, G: nx.Graph, u0: int, v0: int, e0s: list, e1s: list
    ) -> bool:
        """
        Checks to see if the randomly chosen edge is suitable to be swapped. We
        assume that vertices e0[0] -> e1[0] and hence new edges will be created
        between u0 and all e1[1] neighbours; and, between v0 and e0[1] neighbours.
        Edges are not suitable if they don't have the same topology, do have the same
        motif id or if the target edge is already in the graph (possibly in another motif).

        :param G: nx.Graph
        :param u0: focal vertex on left
        :param v0: focal vertex on right
        :param e0s: edges involving focal vertex on left, ebunch
        :param e1s: edges involving focal vertex on right, ebunch
        :return bool: True if edges are OK to be swapped
        """
        if len(e0s) != len(e1s):
            # probably a mismatch during GCM algorithm, selfloops, multi-edges ...
            self._logger.debug(
                f"MarkovChainMonteCarlo - paired corners not equal size {len(e0s)} and {len(e1s)},\
                    {e0s} and {e1s}."
            )
            return False

        hashmap_e0s: dict = self.get_hashmap(G, e0s)
        hashmap_e1s: dict = self.get_hashmap(G, e1s)

        if hashmap_e0s.keys() != hashmap_e1s.keys():
            # check paired corners have same edge-topologies
            self._logger.debug(
                f"MarkovChainMonteCarlo - paired corners not equal topologies {hashmap_e0s} \
                    and {hashmap_e1s}"
            )
            return False

        for key in hashmap_e0s:
            if len(hashmap_e0s[key]) != len(hashmap_e1s[key]):
                # check paired corners have same number of each edge topology, residue
                # of GCM construction
                self._logger.debug(
                    f"MarkovChainMonteCarlo - paired corners not equal number of topologies {hashmap_e0s} \
                        and {hashmap_e1s}"
                )
                return False

        for e0, e1 in zip(e0s, e1s):
            if (
                G.edges[e0][NetworkNames.MOTIF_IDS]
                == G.edges[e1][NetworkNames.MOTIF_IDS]
            ):
                # other chosen vertex is in the same motif
                self._logger.debug(
                    "MarkovChainMonteCarlo - paired corners belong to same motif"
                )
                return False

        # check none of the potential target edges are already present
        # this checks all possible edges that *could* be made.
        for e0 in e0s:
            u1: int = self.get_other_vertex(u0, e0)
            topology: str = G.edges[e0][NetworkNames.TOPOLOGY]
            lst: list = hashmap_e1s[topology]
            for e1 in lst:
                v1: int = self.get_other_vertex(v0, e1)
                if G.has_edge(u0, v1) or G.has_edge(v0, u1):
                    self._logger.debug(
                        "MarkovChainMonteCarlo - target edges already in network"
                    )
                    return False
        return True

    def get_other_vertex(self, u: int, e: tuple) -> int:
        if e[0] == u:
            return e[1]
        elif e[1] == u:
            return e[0]
        else:
            raise ErrorMarkovChainMonteCarloRewiring(
                f"Error: MarkovChainMonteCarlo - get_other_node() vertex {u} not \
                in edge {e}"
            )

    def get_joint_excess_degree_key(self, G: nx.Graph, e: tuple, index: int) -> tuple:
        """
        Returns the joint excess degree key tuple for edge e in network G.
        :param G: nx.Graph
        :param e: edge
        :param index: index in joint degree of edge topology
        :returns key: joint excess joint degree key tuple
        """
        us_jd: list = [list(G.nodes[u][NetworkNames.JOINT_DEGREE]) for u in e]

        us_excess_jd: list = []
        for jd in us_jd:
            jd[index] -= 1
            us_excess_jd.append(tuple(jd))

        return us_excess_jd[0] + us_excess_jd[1]

    def get_swapped_joint_excess_degree_key(
        self, G: nx.Graph, e0: tuple, e1: tuple, u0: int, v0: int, index: int
    ) -> JointExcessJointDegreeKeysView:
        """
        Returns the swapped joint excess degree key tuple for two edges e0, e1 in network G.
        :param G: nx.Graph
        :param e0: edge
        :param e0: edge
        :param u0: left vertex
        :param v0: right vertex
        :param index: index in joint degree of edge topology
        :returns JointExcessJointDegreeKeysView: object to view excess joint degree tuples
        """
        u1: int = self.get_other_vertex(u0, e0)
        v1: int = self.get_other_vertex(v0, e1)

        # pull the joint degrees of the 4 vertices
        us: list = [u0, u1, v0, v1]
        us_jd: list = [list(G.nodes[u][NetworkNames.JOINT_DEGREE]) for u in us]

        # grab their excess degrees in topology `index`
        us_excess_jd: list = []
        for jd in us_jd:
            jd[index] -= 1
            us_excess_jd.append(tuple(jd))

        return JointExcessJointDegreeKeysView(us_excess_jd)

    def append_proposal_edges(
        self, G: nx.Graph, u0: int, old_edge: tuple, new_edge: tuple
    ) -> None:
        """
        Appends edge `e` (tuple) to propsal edges for this trial. Retains the
        focal vertex as first element (u0, v1).
        :param G nx.Graph:
        :param old_edge tuple: edge
        :param new_edge tuple: edge
        """
        p = ProposalEdge()
        p._topology = G.edges[old_edge][NetworkNames.TOPOLOGY]
        p._motif_id = G.edges[old_edge][NetworkNames.MOTIF_IDS]
        v1: int = self.get_other_vertex(u0, new_edge)
        p._new_edge = (u0, v1)
        self._proposal_edges.append(p)

    @MarkovChainMonteCarlo.proposal_efficiency
    def swap_condition(
        self, G: nx.Graph, e0s: list, e1s: list, u0: int, v0: int
    ) -> bool:
        """
        Evaluates the Metropolis condition for the algorithm.

        .. math::
            \\pi=prod_{\\nu}\\prod_i \\frac{e_{tau,k_{u_1},\\nu,k_{v_i}}e_{\tau,k_{v_1},
            nu,k_{v_i}}}{e_{\tau,k_{u_1},\nu,k_{u_i}}e_{\tau,k_{v_1},\nu,k_{v_i}}}

        Ensures that corresponding edges are swapped by creating a topology hashmap
        of (v0,v1) edges so that subgraphs with different orbits swap correct edges.
        Once edges are paired during this step, they are saved in
        `self._proposal_edges`.

        :param G nx.Graph: network
        :param e0s list: list of left vertex edge's in motif, ebunch
        :param e1s list: list of right vertex' edge's in motif, ebunch
        :param u0: left vertex
        :param v0: right vertex
        """
        # create hashmap of (v0,v1) edges keyed by their topologies
        top: float = 1.0
        hashmap_e1s: dict = self.get_hashmap(G, e1s)

        # refresh the list of new edges for this trial
        self._proposal_edges = []

        # iterate (u0,u1) edges
        for e0 in e0s:
            topology: str = G.edges[e0][NetworkNames.TOPOLOGY]
            index: int = self._ejks.get_topology_index(topology)
            # pull partner with same topology from hashmap of e1s
            lst: list = hashmap_e1s[topology]

            try:
                # this is where the other edge is selected
                e1 = lst.pop()
            except IndexError as e:
                raise ErrorMarkovChainMonteCarloRewiring(
                    f"IndexError during swap condition {e}"
                )

            # cache out the new edges that this attempt is proposing to add, (u0, v1)
            # and (v0, u1), to self._proposal_edges
            self.append_proposal_edges(G, u0, e0, (u0, self.get_other_vertex(v0, e1)))
            self.append_proposal_edges(G, v0, e1, (v0, self.get_other_vertex(u0, e0)))

            # create an excess degree key view object
            key_view: JointExcessJointDegreeKeysView = (
                self.get_swapped_joint_excess_degree_key(G, e0, e1, u0, v0, index)
            )

            u0v1_key: tuple = key_view.get_u0v1()
            v0u1_key: tuple = key_view.get_v0u1()
            u0u1: tuple = key_view.get_u0u1()
            u1u0: tuple = key_view.get_u1u0()
            v0v1: tuple = key_view.get_v0v1()
            v1v0: tuple = key_view.get_v1v0()

            # test if swapped keys change anything, if not, silent return False.
            starting_keys: set = {u0u1, u1u0, v0v1, v1v0}
            if set([u0v1_key, v0u1_key]).issubset(starting_keys):
                return False

            try:
                top *= (
                    self._ejks.ejks[topology][u0v1_key]
                    * self._ejks.ejks[topology][v0u1_key]
                )
            except KeyError as e:
                # guard due to GCM algorithm honouring the handshaking lemma by inserting
                # additional motif. Ensure this is not statistically significant, unlikely
                # to fail here due to the previous checks ...
                self._logger.debug(
                    f"MarkovChainMonteCarlo - KeyError during swap condition numerator: {e}"
                )
                return False

            if top == 0.0:
                return False

        bottom: float = 1.0
        # unlike the numerator, the order of iteration of edges in e0s, e1s is unimportant
        # since we are multiplying the all together anyway. However, we should pull the
        # topology for each edge seperately.
        for e0, e1 in zip(e0s, e1s):

            left_topology: str = G.edges[e0][NetworkNames.TOPOLOGY]
            index: int = self._ejks.get_topology_index(left_topology)
            key_e0 = self.get_joint_excess_degree_key(G, e0, index)

            right_topology: str = G.edges[e1][NetworkNames.TOPOLOGY]
            index: int = self._ejks.get_topology_index(right_topology)
            key_e1 = self.get_joint_excess_degree_key(G, e1, index)

            try:
                bottom *= (
                    self._ejks.ejks[left_topology][key_e0]
                    * self._ejks.ejks[right_topology][key_e1]
                )
            except KeyError as e:
                # guard due to GCM algorithm honouring the handshaking lemma by inserting
                # additional motif. Ensure this is not statistically significant.
                self._logger.debug(
                    f"MarkovChainMonteCarlo - KeyError during swap condition denominator: {e}"
                )
                return False

        if bottom == 0.0:
            raise ErrorMarkovChainMonteCarloRewiring("swap_condition() divide by zero")

        value: float = (top + 0.0) / bottom
        if value > random.random():
            return True

        return False

    def rewire(self) -> nx.Graph:
        """
        Entrypoint to MCMC rewiring algorithm
        """
        # pull a copy of the network, careful still has edge data?
        G: nx.Graph = self._network.G.copy()

        number_of_edges: int = G.number_of_edges()

        # create a set based on AVL tree. Ensure to order
        # the edge to avoid duplicate entries for a given edge
        EdgeSet = DrawSet()
        for e in G.edges():
            EdgeSet.add(tuple(sorted(e)))

        convergence_count: int = 0
        while convergence_count <= self._convergence_limit:

            if convergence_count % 50 == 0 and self._proposal_count != 0:
                self._acceptance_ratio.append(
                    float(self._proposals_accepted) / float(self._proposal_count)
                )

            # pick an edge at random.
            e0: tuple = EdgeSet.draw()
            u0: int = e0[0]
            u_edges_in_motif: list = self.get_all_edges(G, u0, e0)

            search_count: int = 0
            while search_count <= self._search_limit:

                # choose another edge at random.
                e1: tuple = EdgeSet.draw()

                if (
                    G.edges[e1][NetworkNames.TOPOLOGY]
                    != G.edges[e0][NetworkNames.TOPOLOGY]
                ):
                    continue

                v0: int = e1[0]
                v_edges_in_motif: list = self.get_all_edges(G, v0, e1)
                if self.is_edge_choice_suitable(
                    G, u0, v0, u_edges_in_motif, v_edges_in_motif
                ):
                    break

                search_count += 1

            if search_count >= self._search_limit:
                self._logger.debug(
                    f"Error: MarkovChainMonteCarlo - rewire() could not find another edge to swap \
                    with {G.edges[e0]}"
                )
                continue

            # we now have two candidate corners of two seperate motifs of the same topology that
            # can be swapped next, evaluate the swap condition and then swap the edges.
            if self.swap_condition(G, u_edges_in_motif, v_edges_in_motif, u0, v0):

                convergence_count += 1

                # add the new proposal edges for both sides
                for pe in self._proposal_edges:
                    e: tuple = pe._new_edge
                    if G.has_edge(*e):
                        raise ErrorMarkovChainMonteCarloRewiring(
                            f"Edge {e} already present in network!"
                        )
                    G.add_edge(*e)
                    G.edges[e][NetworkNames.TOPOLOGY] = pe._topology
                    G.edges[e][NetworkNames.MOTIF_IDS] = pe._motif_id
                    EdgeSet.add(tuple(sorted(e)))

                # remove old edges
                for e0, e1 in zip(u_edges_in_motif, v_edges_in_motif):
                    G.remove_edge(*e0)
                    G.remove_edge(*e1)

                    # EdgeSet expects ordered tuples
                    EdgeSet.remove(tuple(sorted(e0)))
                    EdgeSet.remove(tuple(sorted(e1)))

                if G.number_of_edges() != number_of_edges:
                    raise ErrorMarkovChainMonteCarloRewiring(
                        f"Error: MarkovChainMonteCarlo - rewire() edge count not preserved\
                        by swap {number_of_edges} vs {G.number_of_edges()}"
                    )

        return G

    @property
    def network(self) -> Network:
        return self._network

    @network.setter
    def network(self, value: Network) -> None:
        self._network = value

    @property
    def convergence_limit(self) -> int:
        return self._convergence_limit

    @convergence_limit.setter
    def convergence_limit(self, value: int) -> None:
        self._convergence_limit = value

    @property
    def search_limit(self) -> int:
        return self._search_limit

    @search_limit.setter
    def search_limit(self, value: int) -> None:
        self._search_limit = value

    @property
    def ejks(self) -> JointExcessJointDegreeMatrices:
        return self._ejks

    @ejks.setter
    def ejks(self, value: JointExcessJointDegreeMatrices) -> None:
        self._ejks = value
