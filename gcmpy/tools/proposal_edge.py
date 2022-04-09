class ProposalEdge:
    def __init__(self):
        self._topology: str = None
        self._motif_id: int = None
        self._new_edge: tuple = None

    @property
    def topology(self) -> str:
        return self._topology

    @topology.setter
    def topology(self, value: str) -> None:
        self._topology = value

    @property
    def motif_id(self) -> int:
        return self._motif_id

    @motif_id.setter
    def motif_id(self, value: int) -> None:
        self._motif_id = value

    @property
    def new_edge(self) -> tuple:
        return self._new_edge

    @new_edge.setter
    def new_edge(self, value: tuple) -> None:
        self._new_edge = value
