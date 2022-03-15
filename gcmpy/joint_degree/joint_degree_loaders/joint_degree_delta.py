
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_split_degree import JointDegreeSplitDegree
from gcmpy.joint_degree.joint_degree_types import JointDegreeType

class JointDegreeDelta(JointDegreeSplitDegree):
    
    _type: str = JointDegreeType.DELTA

    def __init__(self, params: dict):
        self._target_k = params["target_k"]
        self._fp = params['fp']
        self._probs = params["probs"]
        self._motif_sizes = params["motif_sizes"]
        if "low_high_degree_bound" in params:
            self._low_high_degree_bound = params["low_high_degree_bound"]

        self.create_jdd()
        
    def create_jdd(self) -> None:
        self._jdd = {}
        for k in range(
            self._low_high_degree_bound[0], self._low_high_degree_bound[1]
            ):
            zeros = [0]*len(self._motif_sizes)
            if k != self._target_k:
                zeros[0] = k
                self._jdd[tuple(zeros)] = self._fp(k)
            else:
                self.resolve_degree(k, self._fp(k))
        self.normalise_jdd()