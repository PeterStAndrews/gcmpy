
from epydemic import DrawSet
import random
import networkx as nx
from logging import Logger

from gcmpy.names.tools_names import ToolsNames
from gcmpy.network.network import Network
from gcmpy.names.network_names import NetworkNames
from gcmpy.tools.joint_excess_joint_degree_matrices import JointExcessJointDegreeMatrices


class MarkovChainMonteCarlo:
    """
    Stochastically rewires `_network` so that the target correlations match `ejks`.
    """
    _network: Network = None
    _convergence_limit: int = None
    _search_limit: int = None
    _ejks: JointExcessJointDegreeMatrices = None

    _logger: Logger = Logger("MCMC")

    def __init__(self, params: dict):
        try:
            self._network: Network = params[ToolsNames.NETWORK]
            self._ejks: JointExcessJointDegreeMatrices = params[ToolsNames.EJKS]
            self._convergence_limit: int = params[ToolsNames.CONVERGENCE_LIMIT]
            self._search_limit: int = params[ToolsNames.SEARCH_LIMIT]
        except Exception as e:
            raise(f"Error instantiating {self.__class__.__name__}: {e}")

    def get_all_edges(self, G: nx.Graph, u0: int, edge: tuple) -> list:
        """
        Returns all edges that have a given motif id as list of ordered tuples
        such that e = (u0,e1).
        :param G: nx.Graph
        :param e: edge in motif
        :returns list: all edges in motif involving focal vertex u0
        """
        es = []
        motif_id: int = G.edges[edge][NetworkNames.MOTIF_IDS]
        for e in G.edges(u0):
            if G.edges[e][NetworkNames.MOTIF_IDS] == motif_id:
                if e[0] == u0:
                    es.append(e)
                else:
                    es.append((e[1], e[0]))
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

    def is_edge_choice_suitable(self, G: nx.Graph, u0: int, v0: int, e0s: list, e1s: list) -> bool:
        """
        Checks to see if the randomly chosen edge is suitable to be swapped. We
        assume that vertices e0[0] -> e1[0] and hence new edges will be created
        between u0 and all e1[1] neighbours; and, between v0 and e0[1] neighbours.
        Edges are not suitable if they don't have the same topology, do have the same
        motif id or if the target edge is already in the graph (possibly in another motif).

        :param G: nx.Graph
        :param u0: focal vertex on left
        :param v0: focal vertex on right
        :param e0s: edges involving focal vertex on left
        :param e1s: edges involving focal vertex on right
        :return bool: True if edges are OK to be swapped
        """
        if len(e0s) != len(e1s):
            self._logger.error(
                f"MarkovChainMonteCarlo - paired corners not equal size {len(e0s)} and {len(e1s)},\
                    {e0s} and {e1s}."
            )
            return False

        hashmap_e0s: dict = self.get_hashmap(G, e0s)
        hashmap_e1s: dict = self.get_hashmap(G, e1s)

        # check paired corners have same edge-topologies
        if hashmap_e0s.keys() != hashmap_e1s.keys():
            self._logger.error(
                f"MarkovChainMonteCarlo - paired corners not equal topologies {hashmap_e0s} \
                    and {hashmap_e1s}"
            )
            return False

        # check paired corners have same number of each edge topology
        for key in hashmap_e0s:
            if len(hashmap_e0s[key]) != len(hashmap_e1s[key]):
                self._logger.error(
                    f"MarkovChainMonteCarlo - paired corners not equal number of topologies {hashmap_e0s} \
                        and {hashmap_e1s}"
                )
                return False

        for e0, e1 in zip(e0s, e1s):
            if G.edges[e0][NetworkNames.MOTIF_IDS] == G.edges[e1][NetworkNames.MOTIF_IDS]:
                self._logger.debug(
                    "MarkovChainMonteCarlo - paired corners belong to same motif"
                )
                return False
            elif G.has_edge(u0, e1[1]) or G.has_edge(v0, e0[1]):
                self._logger.debug(
                    "MarkovChainMonteCarlo - target edges already in network"
                )
                return False
        return True

    def get_other_node(self, u: int, e: tuple) -> int:
        if e[0] == u:
            return e[1]
        elif e[1] == u:
            return e[0]
        else:
            raise f"Error: MarkovChainMonteCarlo - get_other_node() vertex {u} not in edge {e}"

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
            self, G: nx.Graph, e0: tuple, e1: tuple, u0: int, v0: int, index: int) -> tuple:
        """
        Returns the swapped joint excess degree key tuple for two edges e0, e1 in network G. 
        :param G: nx.Graph
        :param e0: edge
        :param e0: edge
        :param u0: left vertex
        :param v0: right vertex
        :param index: index in joint degree of edge topology
        :returns key pair: joint excess joint degree key tuples for (u0,v1) and (v0,u1)
        """
        u1: int = self.get_other_node(u0, e0)
        v1: int = self.get_other_node(v0, e1)
        
        # pull the joint degrees of the 4 vertices
        us: list = [u0,u1,v0,v1]
        us_jd: list = [list(G.nodes[u][NetworkNames.JOINT_DEGREE])for u in us]

        # grab their excess degrees in topology `index`
        us_excess_jd: list = []
        for jd in us_jd:
            jd[index] -= 1
            us_excess_jd.append(tuple(jd))

        return us_excess_jd[0] + us_excess_jd[3], us_excess_jd[1] + us_excess_jd[2] 

    def swap_condition(self, G: nx.Graph, e0s: list, e1s: list, u0: int, v0: int) -> bool:
        """
        Evaluates the Metropolis condition for the algorithm. 

        .. math::
            \pi=\prod_{\nu}\prod_i\frac{e_{\tau,k_{u_1},\nu,k_{v_i}}e_{\tau,k_{v_1},
            \nu,k_{v_i}}}{e_{\tau,k_{u_1},\nu,k_{u_i}}e_{\tau,k_{v_1},\nu,k_{v_i}}}

        Ensures that corresponding edges are swapped by creating a topology hashmap
        of (v0,v1) edges so that subgraphs with different orbits swap correct edges.

        :param G nx.Graph: network
        :param e0s list: list of left vertex edge's in motif
        :param e1s list: list of right vertex' edge's in motif
        :param u0: left vertex
        :param v0: right vertex
        """
        # create hashmap of (v0,v1) edges keyed by their topologies
        top: float = 1.0
        hashmap_e1s: dict = self.get_hashmap(G, e1s)

        # iterate (u0,u1) edges
        for e0 in e0s:
            topology: str = G.edges[e0][NetworkNames.TOPOLOGY]
            index: int = self._ejks.get_topology_index(topology)
            # pull partner with same topology from hashmap of e1s
            lst: list = hashmap_e1s[topology]
            e1 = lst.pop()
            key1, key2 = self.get_swapped_joint_excess_degree_key(
                G, e0, e1, u0, v0, index
            )
            try:
                # guard due to GCM algorithm honouring the handshaking lemma by inserting
                # additional motif. Ensure this is not statistically significant.
                top *= self._ejks.ejks[topology][key1] * self._ejks.ejks[topology][key2]
            except KeyError as e:
                self._logger.error(
                    f"MarkovChainMonteCarlo - KeyError during swap condition numerator: {e}"
                )
                return False

            if top == 0.0:
                return False

        bottom: float = 1.0
        for e0, e1 in zip(e0s, e1s):
            
            topology: str = G.edges[e0][NetworkNames.TOPOLOGY]
            topology_e1: str = G.edges[e1][NetworkNames.TOPOLOGY]

            if topology != topology_e1:
                raise f"Error: MarkovChainMonteCarlo - swap_condition() edge topologies {topology} != {topology_e1}"
            
            index: int = self._ejks.get_topology_index(topology)
            key_e0 = self.get_joint_excess_degree_key(G, e0, index)
            key_e1 = self.get_joint_excess_degree_key(G, e1, index)

            try:
                # guard due to GCM algorithm honouring the handshaking lemma by inserting
                # additional motif. Ensure this is not statistically significant.
                bottom *= self._ejks.ejks[topology][key_e0] * self._ejks.ejks[topology][key_e1]
            except KeyError as e:
                self._logger.error(
                    f"MarkovChainMonteCarlo - KeyError during swap condition denominator: {e}"
                )
                return False

        if bottom == 0.0:
            raise "Error: MarkovChainMonteCarlo - swap_condition() divide by zero"

        if (top + 0.0) / bottom > random.random():
            return True

        return False

    def rewire(self) -> None:

        # pull a copy of the network, careful still has edge data?
        G: nx.Graph = self._network._G

        number_of_edges: int = G.number_of_edges()

        # create a set based on AVL tree. Ensure to order 
        # the edge to avoid duplicate entries for a given edge
        EdgeSet = DrawSet()
        for e in G.edges():
            if e[0] <= e[1]:
                EdgeSet.add(e)
            else:
                EdgeSet.add((e[1], e[0]))

        convergence_count: int = 0
        while convergence_count <= self._convergence_limit:

            # pick an edge at random.
            e0: tuple = EdgeSet.draw()
            u0: int = e0[0]
            u_edges_in_motif: list = self.get_all_edges(G, u0, e0)

            search_count: int = 0
            while search_count <= self._search_limit:

                # choose another edge at random.
                e1: tuple = EdgeSet.draw()

                if G.edges[e1][NetworkNames.TOPOLOGY] != G.edges[e0][NetworkNames.TOPOLOGY]:
                    continue

                v0: int = e1[0]
                v_edges_in_motif: list = self.get_all_edges(G, v0, e1)
                if self.is_edge_choice_suitable(G, u0, v0, u_edges_in_motif, v_edges_in_motif):
                    break

                search_count += 1

            if search_count >= self._search_limit:
                self._logger.error(
                    f"Error: MarkovChainMonteCarlo - rewire() could not find another edge to swap with {G.edges[e0]}"
                )
                continue

            # we now have two candidate corners of two seperate motifs that can be swapped
            # next, evaluate the swap condition and then swap the edges.
            if self.swap_condition(G, u_edges_in_motif, v_edges_in_motif, u0, v0):

                convergence_count += 1

                # add the new edges between v0 & u1 vertices retaining topology of old (u0,u1) edge
                for e in u_edges_in_motif:
                    topology: str = G.edges[e][NetworkNames.TOPOLOGY]
                    motif_id: int = G.edges[e][NetworkNames.MOTIF_IDS]
                    u1: int = self.get_other_node(u0, e)
                    G.add_edge(v0, u1)
                    G.edges[(v0, u1)][NetworkNames.TOPOLOGY] = topology
                    G.edges[(v0, u1)][NetworkNames.MOTIF_IDS] = motif_id
                    if v0 <= u1:
                        EdgeSet.add((v0,u1))
                    else:
                        EdgeSet.add((u1,v0))

                for e in v_edges_in_motif:
                    topology: str = G.edges[e][NetworkNames.TOPOLOGY]
                    motif_id: int = G.edges[e][NetworkNames.MOTIF_IDS]
                    v1: int = self.get_other_node(v0, e)
                    G.add_edge(u0, v1)
                    G.edges[(u0, v1)][NetworkNames.TOPOLOGY] = topology
                    G.edges[(u0, v1)][NetworkNames.MOTIF_IDS] = motif_id
                    if u0 <= v1:
                        EdgeSet.add((u0,v1))
                    else:
                        EdgeSet.add((v1,u0))

                # remove old edges
                for e0, e1 in zip(u_edges_in_motif, v_edges_in_motif):
                    G.remove_edge(*e0)
                    G.remove_edge(*e1)

                    # EdgeSet expects ordered tuples
                    if e0[0] > e0[1]:
                        EdgeSet.remove((e0[1],e0[0]))
                    else:
                        EdgeSet.remove(e0)
                    
                    if e1[0] > e1[1]:
                        EdgeSet.remove((e1[1],e1[0]))
                    else:
                        EdgeSet.remove(e1)
                    

                if G.number_of_edges() != number_of_edges:
                    raise f"Error: MarkovChainMonteCarlo - rewire() edge count not preserved\
                         by swap {number_of_edges} vs {G.number_of_edges()}"
