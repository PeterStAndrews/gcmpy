from abc import ABC, abstractmethod
from typing import Any, Generator

from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames


class GCMAlgorithm(ABC):
    """
    Abstract base class for GCM algorithm. Expects a parameter dict with keys:
    :param motif_sizes: list of ints that indicate the number of vertices in each motif
    :param build_functions: callbacks that accept list of vertices and return edges
    :param edge_names: list of names for edge topologies
    """

    def __init__(self, params: dict):
        self._motif_sizes: list = []  # list of number of vertices in each motif
        self._build_functions: list = []  # list of callbacks for motif construction
        self._edge_names: list = []  # list of edge topology names

        try:
            self._motif_sizes = params[GCMAlgorithmNames.MOTIF_SIZES]
            self._build_functions = params[GCMAlgorithmNames.BUILD_FUNCTIONS]
            self._edge_names = params[GCMAlgorithmNames.EDGE_NAMES]
        except Exception as e:
            raise (f"Error in {self.__class__.__name__}: {e}")

    def __new__(cls, *args, **kwargs):
        if cls is GCMAlgorithm:
            raise TypeError(
                "The GCMAlgorithm class is abstract and may not be instantiated"
            )
        return object.__new__(cls)

    @abstractmethod
    def random_clustered_graph(self, jds: list) -> Any:
        raise NotImplementedError(
            "Error attempting to call virtual method on GCMAlgorithm: random_clustered_graph"
        )

    def infinite_sequence(self) -> Generator:
        num = 0
        while True:
            yield num
            num += 1
