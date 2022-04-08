from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_type import JointDegreeType
from gcmpy.names.joint_degree_names import JointDegreeNames


class JointDegreeCover(JointDegree):

    _type: str = JointDegreeType.COVER

    def __init__(self, params: dict):
        try:
            self._cover = params[JointDegreeNames.COVER]
        except Exception as e:
            raise (f"Error instantiating {self.__class__.__name__}: {e}")

        self._motif_sizes = sorted(list(set([len(c) for c in self._cover])))
        self.create_jdd()

    def create_jdd(self) -> None:

        vertex_ids = list(set([vertex for clique in self._cover for vertex in clique]))

        zero_index = 0
        if min(vertex_ids) != zero_index:
            zero_index = 1

        largest_clique = len(max(self._cover, key=len))

        jds = []
        for _ in range(len(vertex_ids)):
            jd = [0] * largest_clique
            jds.append(jd)

        for c in self._cover:
            clique_size = len(c)
            for vertex in c:
                jds[vertex - zero_index][clique_size - 1] += 1

        # iterate each column of the jds and record the index if all zeros
        indxs = [i for i, top in enumerate(zip(*jds)) if not any(top)]

        # use the indexes of the zero columns to remove
        for i in indxs:
            for jd in jds:
                del jd[i]

        # convert jds to jdd
        self.convert_jds_to_jdd(jds)

    @property
    def cover(self) -> list:
        return self._cover

    @cover.setter
    def cover(self, value: list) -> None:
        self._cover = value
