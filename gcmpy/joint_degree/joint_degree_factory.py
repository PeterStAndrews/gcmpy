from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_type import JointDegreeType
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import (
    JointDegreeManual,
)
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_empirical import (
    JointDegreeEmpirical,
)
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_function import (
    JointDegreeFunction,
)
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_marginal import (
    JointDegreeMarginal,
)
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_split_degree import (
    JointDegreeSplitDegree,
)
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_delta import JointDegreeDelta
from gcmpy.joint_degree.joint_degree_loaders.joint_degree_cover import JointDegreeCover


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
