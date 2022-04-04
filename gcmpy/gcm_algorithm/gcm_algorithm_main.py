from gcmpy.gcm_algorithm.gcm_algorithm_types import GCMAlgorithmTypes
from gcmpy.gcm_algorithm.gcm_algorithm_factory import GCMAlgorithmFactory
from gcmpy.gcm_algorithm.gcm_algorithm import GCMAlgorithm
from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames


class GCMAlgorithmMain:
    """
    Determines the type of GCM algorithm subclass to return. Expects `params`
    dict to contain key `GCM_type` among others additionally required
    by the subclasses.
    :method load_gcm_algorithm: returns the subclass of ABC `GCMAlgorithm`.
    This will raise an error if the params dict is not correctly set up.
    """

    @staticmethod
    def load_gcm_algorithm(params: dict) -> GCMAlgorithm:
        try:
            input_type = GCMAlgorithmTypes(params[GCMAlgorithmNames.GCM_TYPE])
            loader: GCMAlgorithm = GCMAlgorithmFactory.resolve_algorithm(
                input_type, params
            )
            return loader
        except Exception as e:
            raise (f"Error instantiating GeneralisedConfigurationModel: {e}")
