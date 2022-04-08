from itertools import product
import numpy as np
import random

from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_type import JointDegreeType
from gcmpy.names.joint_degree_names import JointDegreeNames


class JointDegreeMarginal(JointDegree):

    """
    Merge uncorrelated marginals in each topology from analytical data together to create
    self._jdd. If using a direct method, all possible joint degree tuples are evaluated;
    however, for large varience in the allowed degrees this method is slow. Instead, we can
    choose to sample the analytical functions by setting use_sampling which draws n_samples
    weighted samples from each marginal function and pieces them together.
    :param arr_fp: array of callbacks
    :param motif_sizes: list of ints for number of vertices in each motif
    :param hi_lo_degree_bounds: list of tuples (int,int) for kmin,kmax per topology
    :param use_sampling: bool to use sampling or direct approach
    :param n_samples: number of samples if not direct
    """

    _type: str = JointDegreeType.MARGINAL

    def __init__(self, params: dict):
        self._arr_fp: list[callable] = []
        self._low_high_degree_bounds: tuple = None
        self._use_sampling: bool = False
        self._n_samples: int = 100000

        try:
            self._motif_sizes = params[JointDegreeNames.MOTIF_SIZES]
            self._arr_fp = params[JointDegreeNames.ARR_FP]
            self._low_high_degree_bounds = params[
                JointDegreeNames.LOW_HIGH_DEGREE_BOUND
            ]
        except Exception as e:
            raise (f"Error instantiating {self.__class__.__name__}: {e}")

        # default not to use sampling method
        if JointDegreeNames.USE_SAMPLING in params:
            self._use_sampling = params[JointDegreeNames.USE_SAMPLING]
        if JointDegreeNames.N_SAMPLES in params:
            self._n_samples = params[JointDegreeNames.N_SAMPLES]
        self.create_jdd()

    def create_jdd(self) -> None:
        if not self._use_sampling:
            self.create_jdd_directly()
        else:
            self.create_jdd_by_sampling()

    def create_jdd_directly(self) -> None:
        """
        Create self._jdd by enumerating probability of all possible joint degrees from
        analytical functions.
        """
        self._jdd = dict((key, 0.0) for key in self.generate_all_joint_degrees())
        for key in self._jdd:
            self._jdd[key] = self.evaluate_prob_of_joint_degree(key)
        self.normalise_jdd()

    def generate_all_joint_degrees(self) -> list:
        """
        Generate all possible joint degrees from a range of min/max degree of each motif.
        :return jd: list of joint degrees
        """
        ks = []
        for kmin, kmax in self._low_high_degree_bounds:
            ks.append([k for k in range(kmin, kmax)])
        return list(product(*ks))

    def evaluate_prob_of_joint_degree(self, joint_degree: list) -> float:
        """
        Evaluate the joint probability of a joint degree using function callbacks.
        :param joint_degree: list of joint degreees
        :returns prod: probability of joint degree
        """
        prod: float = 1.0
        for i, deg in enumerate(joint_degree):
            prod *= self._arr_fp[i](deg)
        return prod

    def draw_from_analytical_joint(self) -> list:
        """Draw `self._n_samples' from a marginal distribution `self._arr_fp' along each dimension.
        :returns jds: samples from marginal distributions"""
        ret = []
        for i in range(len(self._low_high_degree_bounds)):
            kmin, kmax = self._low_high_degree_bounds[i]
            ks = [k for k in range(kmin, kmax + 1)]  # possible degrees
            pks = [self._arr_fp[i](k) for k in ks]  # degree weights
            ret.append(
                random.choices(ks, pks, k=self._n_samples)
            )  # sample this dimension
        return [
            tuple(jd) for jd in np.column_stack(ret).tolist()
        ]  # return sampled degrees

    def create_jdd_by_sampling(self) -> None:
        """
        Draws jds samples from analytical marginal functions and converts to
        a joint degree distribution.
        """
        self.convert_jds_to_jdd(self.draw_from_analytical_joint())
