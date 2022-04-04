from gcmpy.gcm_algorithm.gcm_algorithm_types import GCMAlgorithmTypes
from gcmpy.gcm_algorithm.gcm_algorithm import GCMAlgorithm

from gcmpy.gcm_algorithm.gcm_algorithm_fast import GCMAlgorithmFast
from gcmpy.gcm_algorithm.gcm_algorithm_network import GCMAlgorithmNetwork
from gcmpy.gcm_algorithm.gcm_algorithm_custom_motifs import GCMAlgorithmCustomMotifs


class GCMAlgorithmFactory:
    @staticmethod
    def resolve_algorithm(type: GCMAlgorithmTypes, params: dict) -> GCMAlgorithm:
        if type == GCMAlgorithmTypes.FAST:
            return GCMAlgorithmFast(params)
        elif type == GCMAlgorithmTypes.NETWORK:
            return GCMAlgorithmNetwork(params)
        elif type == GCMAlgorithmTypes.MOTIFS:
            return GCMAlgorithmCustomMotifs(params)
        else:
            raise ("Error: unknown algorithm in GCMAlgorithmFactory: resolve_algorithm")
