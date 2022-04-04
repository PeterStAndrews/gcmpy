from enum import Enum


class JointDegreeType(Enum):
    MANUAL = "manual"
    EMPIRICAL = "empirical"
    JOINT_FUNCTION = "function"
    MARGINAL = "marginal"
    SPLIT_DEGREE = "split_degree"
    DELTA = "delta"
    COVER = "cover"
    UNDEFINED = "undefined"
