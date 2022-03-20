
from abc import ABC, abstractmethod
from typing import Any

class GCMAlgorithm(ABC):
    """
    Abstract base class for GCM algorithm. Expects a parameter dict with keys:
    :param motif_sizes: list of ints that indicate the number of nodes in each motif
    :param build_functions: callbacks that accept list of nodes and return edges
    :param edge_names: list of names for edge topologies
    """
    _motif_sizes: list                      # list of number of nodes in each motif
    _build_functions: list                  # list of callbacks for motif construction 
    _edge_names: list                       # list of edge topology names
    
    def __init__(self, params: dict):
        try:
            self._motif_sizes = params["motif_sizes"]
            self._build_functions = params["build_functions"]
            self._edge_names = params["edge_names"]
        except Exception as e:
            raise (f"Error in GCMAlgorithm: {e}")

    def __new__(cls,*args,**kwargs):
        if cls is GCMAlgorithm:
            raise TypeError(
                "The GCMAlgorithm class is abstract and may not be instantiated"
            )
        return object.__new__(cls)

    @abstractmethod
    def random_clustered_graph(self, jds: list) -> Any:
        raise NotImplementedError (
            "Error attempting to call virtual method on GCMAlgorithm: random_clustered_graph"
        )