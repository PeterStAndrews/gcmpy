from gcmpy.names.tools_names import ToolsNames


class JointExcessJointDegreeMatrices:
    def __init__(self, params=None):
        self._ejks: dict[dict] = {}  # dict of ejk dicts/matrices keyed by edge topology
        self._excess_degree_keys: dict[
            list
        ] = {}  # dict of lst of excess deg-tuples keyed by edge topology
        self._topology_names: list[str] = []  # list of str of topology names, ordered

        if params is not None:
            self._ejks = params[ToolsNames.EJKS]
            self._topology_names = params[ToolsNames.EDGE_NAMES]
            self.get_excess_degree_keys()

    def get_topology_index(self, topology: str) -> int:
        for i, name in enumerate(self._topology_names):
            if name == topology:
                return i
        raise f"Error: JointExcessJointDegreeMatrices - no topology {name} in \
            names {self._topology_names}"

    def get_excess_degree_keys(self) -> None:
        """
        Method to get excess degree keys from ejk matrices for each topology.
        """
        self._excess_degree_keys: dict[list] = {}
        for topology in self._ejks:

            excess_key_set = set()
            for excess_key in self._ejks[topology]:

                A = list(excess_key)
                B = A[: len(A) // 2]
                C = A[len(A) // 2 :]
                excess_key_set.add(tuple(B))
                excess_key_set.add(tuple(C))

            self._excess_degree_keys[topology] = list(excess_key_set)

    @property
    def excess_degree_keys(self) -> dict:
        return self._excess_degree_keys

    @excess_degree_keys.setter
    def excess_degree_keys(self, value: dict) -> None:
        self._excess_degree_keys = value

    @property
    def ejks(self) -> dict:
        return self._ejks

    @ejks.setter
    def ejks(self, value: dict) -> None:
        self._ejks = value

    @property
    def topology_names(self) -> list:
        return self._topology_names

    @topology_names.setter
    def topology_names(self, value: list) -> None:
        self._topology_names = value
