from gcmpy.joint_degree.joint_degree_loaders.joint_degree_split_degree import (
    JointDegreeSplitDegree,
)
from gcmpy.joint_degree.joint_degree_type import JointDegreeType
from gcmpy.names.joint_degree_names import JointDegreeNames


class JointDegreeDelta(JointDegreeSplitDegree):

    _type: str = JointDegreeType.DELTA

    def __init__(self, params: dict):

        try:
            self._target_k = params[JointDegreeNames.TARGET_K]
            self._fp = params[JointDegreeNames.FP]
            self._probs = params[JointDegreeNames.PROBS]
            self._motif_sizes = params[JointDegreeNames.MOTIF_SIZES]
            self._low_high_degree_bound = params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND]
        except Exception as e:
            raise (f"Error instantiating {self.__class__.__name__}: {e}")

        self.create_jdd()

    def create_jdd(self) -> None:
        self._jdd = {}
        for k in range(self._low_high_degree_bound[0], self._low_high_degree_bound[1]):
            zeros = [0] * len(self._motif_sizes)
            if k != self._target_k:
                zeros[0] = k
                self._jdd[tuple(zeros)] = self._fp(k)
            else:
                self.resolve_degree(k, self._fp(k))
        self.normalise_jdd()
