from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_type import JointDegreeType
from gcmpy.names.joint_degree_names import JointDegreeNames


class JointDegreeSplitDegree(JointDegree):

    _type: str = JointDegreeType.SPLIT_DEGREE

    def __init__(self, params: dict):
        self._fp: callable = None
        self._probs: list = []
        self._low_high_degree_bound: tuple = None

        try:
            self._fp = params[JointDegreeNames.FP]
            self._probs = params[JointDegreeNames.PROBS]
            self._motif_sizes = params[JointDegreeNames.MOTIF_SIZES]
            self._low_high_degree_bound = params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND]
        except Exception as e:
            raise (f"Error instantiating {self.__class__.__name__}: {e}")

        self.create_jdd()

    def create_jdd(self) -> None:
        for k in range(self._low_high_degree_bound[0], self._low_high_degree_bound[1]):
            self.resolve_degree(k, self._fp(k))
        self.normalise_jdd()

    def get_valid_joint_degrees(self, remaining_degree: int, topology: int) -> list:
        """
        Returns a list of tuples by recursion. Only an ordered list of cliques are
        currently supported.
        :param remaining_degree: current free edges that can be partitioned
        :param topology: column index of joint degree tuple
        :returns list of joint degrees
        """
        if topology == 1:
            yield [remaining_degree]
        else:
            for i in range(0, remaining_degree // topology + 1):
                for row in self.get_valid_joint_degrees(
                    remaining_degree - i * topology, topology - 1
                ):
                    yield row + [i]

    def calc_prob_of_joint_degree(self, jd: tuple) -> float:
        """
        Calculates the probability of a joint degree from the input params.
        :param jd: joint degree
        :returns probability: float value
        """
        prod: float = 1.0
        for i, degree in enumerate(jd):
            prod *= pow(self._probs[i], (i + 1) * degree)
        return prod

    def resolve_degree(self, k: int, prob_overall_k: float) -> None:
        """
        creates all valid joint degrees for overall k and their probability
        and updates self._jdd.
        :param k: overall degree
        :param prob_overall_k: float value
        """
        self._jdd = {}

        # get a list of valid joint degrees
        valid_tuples: list = list(self.get_valid_joint_degrees(k, len(self._probs)))

        # calculate the probability of each joint degree tuple
        probabilities = []
        for jd in valid_tuples:
            probabilities.append(self.calc_prob_of_joint_degree(jd))

        # normalise the probabilities to unity
        total = sum(probabilities)
        for i in range(len(probabilities)):
            probabilities[i] /= total

        # add each tuple to self._jdd weighted by probability of overall degree k
        for i, jd in enumerate(valid_tuples):
            self._jdd[tuple(jd)] = prob_overall_k * probabilities[i]
