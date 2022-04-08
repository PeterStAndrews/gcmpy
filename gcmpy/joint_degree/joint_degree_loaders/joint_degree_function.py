from itertools import product

from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_type import JointDegreeType
from gcmpy.names.joint_degree_names import JointDegreeNames


class JointDegreeFunction(JointDegree):
    """
    Multivariate function to evaluate the probability of given joint degree from an analytical
    source. Note, callable self._fp signature must accept a joint degree tuple and return a float.
    :param fp: callback
    :param motif_sizes: list of ints for number of vertices in each motif
    :param hi_lo_degree_bounds: list of tuples (int,int) for kmin,kmax per topology
    :param n_samples: number of samples if not direct
    """

    _type: str = JointDegreeType.JOINT_FUNCTION

    def __init__(self, param: dict):
        self._fp: callable = None
        self._low_high_degree_bounds: tuple = (0, 50)
        try:
            self._motif_sizes = param[JointDegreeNames.MOTIF_SIZES]
            self._fp = param[JointDegreeNames.FP]
            self._low_high_degree_bounds = param[JointDegreeNames.LOW_HIGH_DEGREE_BOUND]
        except Exception as e:
            raise (f"Error instantiating {self.__class__.__name__}: {e}")

        self.create_jdd()

    def create_jdd(self) -> None:
        """
        Evaluates probability directly by generating all possible joint degrees.
        """
        # build list of lists of possible degrees in each dimension
        ks = [
            list(range(kmin, kmax + 1)) for kmin, kmax in self._low_high_degree_bounds
        ]
        # iterate all joint degrees and evaluate the joint degree
        for jd in list(product(*ks)):
            self._jdd[jd] = self._fp(jd)
