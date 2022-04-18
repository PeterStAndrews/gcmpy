from gcmpy import JointDegree
from gcmpy import JointDegreeType
from gcmpy import (
    JointDegreeManual,
)
from gcmpy import (
    JointDegreeEmpirical,
)
from gcmpy import (
    JointDegreeFunction,
)
from gcmpy import (
    JointDegreeMarginal,
)
from gcmpy import (
    JointDegreeSplitDegree,
)
from gcmpy import JointDegreeDelta
from gcmpy import JointDegreeCover


class JointDegreeFactory:
    @staticmethod
    def resolve_joint_degree(type: JointDegreeType, params: dict) -> JointDegree:
        if type == JointDegreeType.MANUAL:
            return JointDegreeManual(params)
        elif type == JointDegreeType.EMPIRICAL:
            return JointDegreeEmpirical(params)
        elif type == JointDegreeType.JOINT_FUNCTION:
            return JointDegreeFunction(params)
        elif type == JointDegreeType.MARGINAL:
            return JointDegreeMarginal(params)
        elif type == JointDegreeType.SPLIT_DEGREE:
            return JointDegreeSplitDegree(params)
        elif type == JointDegreeType.DELTA:
            return JointDegreeDelta(params)
        elif type == JointDegreeType.COVER:
            return JointDegreeCover(params)
        else:
            raise ValueError("Unknown joint degree type")
