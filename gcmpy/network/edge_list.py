class LightWeightEdgeList:
    def __init__(self):
        self._edge_list: list = []
        self._topologies: list = []
        self._joint_degrees: list = []
        self._motif_id: list = []

    @property
    def edge_list(self) -> list:
        return self._edge_list

    @edge_list.setter
    def edge_list(self, value: list) -> None:
        self._edge_list = value

    @property
    def topologies(self) -> list:
        return self._topologies

    @topologies.setter
    def topologies(self, value: list) -> None:
        self._topologies = value

    @property
    def joint_degrees(self) -> list:
        return self._joint_degrees

    @joint_degrees.setter
    def joint_degrees(self, value: list) -> None:
        self._joint_degrees = value

    @property
    def motif_id(self) -> list:
        return self._motif_id

    @motif_id.setter
    def motif_id(self, value: list) -> None:
        self._motif_id = value
