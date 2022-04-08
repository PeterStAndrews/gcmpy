import networkx as nx
from gcmpy.names.network_names import NetworkNames
from gcmpy.names.tools_names import ToolsNames
from gcmpy.tools.joint_excess_joint_degree_matrices import (
    JointExcessJointDegreeMatrices,
)


class JointExcessJointDegree:

    """
    Extracts the joint excess joint degree distribution from a network
    """

    def __init__(self, params: dict):
        self._topology_names: list[str] = None
        self._G: nx.Graph = None
        self._ejks = None
        self._num_edges: dict = {}  # number of edges of i-th topology
        self._degree_keys: list = []  # list of joint degree tuples
        self._excess_degree_keys = {}

        try:
            self._G: nx.Graph = params[ToolsNames.NETWORK]
            self._topology_names: list[str] = params[ToolsNames.EDGE_NAMES]
            self._degree_keys = set(
                [
                    tuple(self._G.nodes[n][NetworkNames.JOINT_DEGREE])
                    for n in self._G.nodes()
                ]
            )
        except Exception as e:
            raise (f"Error in {self.__class__.__name__}: {e}")

        self.resolve_excess_degree_keys()

    def resolve_excess_degree_keys(self) -> None:
        """
        Obtains the excess degree keys from the joint degree keys
        for each topology.
        """
        for i, topology in enumerate(self._topology_names):
            excess_keys_top_i = []
            for jd in self._degree_keys:
                if jd[i] > 0:
                    jd = list(jd)
                    jd[i] -= 1
                    excess_keys_top_i.append(tuple(jd))
            self._excess_degree_keys[topology] = list(set(excess_keys_top_i))

    def count_edge_types(self) -> None:
        """
        Enumerates the number of edges of the i-th topology.
        Assumes that the network edges have attribute `topology`.
        """
        # self._num_edges = {}
        for e in self._G.edges():
            topology: str = self._G.edges[e][NetworkNames.TOPOLOGY]
            self._num_edges[topology] = self._num_edges.get(topology, 0) + 1

    def get_ejk(self, i: int, name: str) -> dict:
        """
        Get mixing patterns for the i-th topology. Routine assumes that
        the network edges have attribute `topology` and that the vertices
        have attribute `joint_degree`.
        :param i int: index of topology i in joint degree
        :param name str: key for edge name
        :returns ejk dict: mixing patterns for this topology
        """
        ejk = {}
        for e in self._G.edges():
            if self._G.edges[e][NetworkNames.TOPOLOGY] == name:

                u, v = e

                u_joint_degree = list(self._G.nodes[u][NetworkNames.JOINT_DEGREE])
                v_joint_degree = list(self._G.nodes[v][NetworkNames.JOINT_DEGREE])

                u_joint_degree[i] -= 1
                u_joint_excess_degree = tuple(u_joint_degree)

                v_joint_degree[i] -= 1
                v_joint_excess_degree = tuple(v_joint_degree)

                key1 = u_joint_excess_degree + v_joint_excess_degree
                key2 = v_joint_excess_degree + u_joint_excess_degree

                if key1 == key2:
                    ejk[key1] = ejk.get(key1, 0) + (1.0 / (self._num_edges[name]))
                else:
                    ejk[key1] = ejk.get(key1, 0) + (1.0 / (self._num_edges[name]))
                    ejk[key2] = ejk.get(key2, 0) + (1.0 / (self._num_edges[name]))

        return ejk

    def get_ejks(self) -> JointExcessJointDegreeMatrices:
        """
        Get mixing patterns for all edge topologies
        :returns self._ejks: JointExcessJointDegreeMatrices.
        """
        self._ejks = JointExcessJointDegreeMatrices()
        self._ejks.excess_degree_keys = self._excess_degree_keys
        self._ejks.topology_names = self._topology_names

        self.count_edge_types()
        for i, topology in enumerate(self._topology_names):
            self._ejks.ejks[topology] = self.get_ejk(i, topology)
        return self._ejks
