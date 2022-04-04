from gcmpy.joint_degree.joint_degree_type import JointDegreeType
from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_factory import JointDegreeFactory
from gcmpy.names.joint_degree_names import JointDegreeNames


class JointDegreeDistribution:
    """
    Determines the type of JointDegree subclass to return. Expects `params`
    dict to contain key `joint_degree_type'.
    :method load_joint_degree: returns a subclass of ABC `JointDegree`. This
    will raise an error if the params dict does not contain the required keys.
    """

    @staticmethod
    def load_joint_degree(params: dict) -> JointDegree:
        try:
            input_type = JointDegreeType(params[JointDegreeNames.JOINT_DEGREE_TYPE])
        except Exception as e:
            raise (f"Error instantiating JointDegreeDistribution: {e}")
        loader: JointDegree = JointDegreeFactory.resolve_joint_degree(
            input_type, params
        )
        loader.create_jdd()
        return loader
