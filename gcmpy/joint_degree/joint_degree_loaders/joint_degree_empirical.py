from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_type import JointDegreeType
from gcmpy.names.joint_degree_names import JointDegreeNames


class JointDegreeEmpirical(JointDegree):

    _type: str = JointDegreeType.EMPIRICAL

    def __init__(self, params: dict):
        self._empirical_jds: list = []
        try:
            self._motif_sizes = params[JointDegreeNames.MOTIF_SIZES]
            self._empirical_jds = params[JointDegreeNames.JDS]
        except Exception as e:
            raise (f"Error instantiating {self.__class__.__name__}: {e}")
        self.create_jdd()

    def create_jdd(self) -> None:
        self.convert_jds_to_jdd(self._empirical_jds)

    @property
    def empirical_jds(self) -> list:
        return self._empirical_jds

    @empirical_jds.setter
    def empirical_jds(self, value: list) -> None:
        self._empirical_jds = value
