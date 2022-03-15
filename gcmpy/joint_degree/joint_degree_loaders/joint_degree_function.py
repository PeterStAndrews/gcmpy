
from itertools import product

from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_types import JointDegreeType

class JointDegreeFunction(JointDegree):
    """
    Multivariate function to evaluate the probability of given joint degree from an analytical source.
    Note, callable self._fp signature must accept a joint degree tuple and return a float.
    :param fp: callback
    :param motif_sizes: list of ints for number of vertices in each motif
    :param hi_lo_degree_bounds: list of tuples (int,int) for kmin,kmax per topology
    :param n_samples: number of samples if not direct
    """
    _type: str = JointDegreeType.JOINT_FUNCTION
    _fp: callable = None 
    _low_high_degree_bounds: tuple = (0,50)

    def __init__(self, param: dict):      
        self._motif_sizes = param["motif_sizes"]
        self._fp = param["fp"]
        # important for powerlaw distribution functions to set high
        if "low_high_degree_bounds" in param:
            self._low_high_degree_bounds = param["low_high_degree_bounds"]
        self.create_jdd()
        
    def create_jdd(self) -> None:
        """
        Evaluates probability directly by generating all possible joint degrees.
        """
        # build list of lists of possible degrees in each dimension
        ks = [list(range(kmin,kmax+1)) for kmin,kmax in self._low_high_degree_bounds]
        # iterate all joint degrees and evaluate the joint degree
        for jd in list(product(*ks)):
            self._jdd[jd] = self._fp(jd)