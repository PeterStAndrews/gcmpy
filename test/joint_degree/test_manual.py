import unittest

from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import (
    JointDegreeManual,
)
from gcmpy.names.joint_degree_names import JointDegreeNames

NETWORK_SIZE: int = 100000


class JDManualTest(unittest.TestCase):
    def test_manual_JDD_single_topology(self):

        # valid input data for manual entry
        params = {}
        params[JointDegreeNames.JDD] = {(1,): 0.2, (2,): 0.5, (3,): 0.1, (5,): 0.2}
        params[JointDegreeNames.MOTIF_SIZES] = [2]

        DegreeDistObj = JointDegreeManual(params)
        n_vertices: int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params[JointDegreeNames.MOTIF_SIZES])):
            self.assertTrue(
                sum([jd[0] for jd in jds]) % params[JointDegreeNames.MOTIF_SIZES][i]
                == 0
            )

    def test_manual_JDD_two_topologies(self):

        params = {}
        params[JointDegreeNames.JDD] = {
            (1, 0): 0.2,
            (2, 1): 0.5,
            (3, 0): 0.1,
            (5, 1): 0.2,
        }
        params[JointDegreeNames.MOTIF_SIZES] = [2, 3]

        DegreeDistObj = JointDegreeManual(params)
        n_vertices: int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params[JointDegreeNames.MOTIF_SIZES])):
            self.assertTrue(
                sum([jd[i] for jd in jds]) % params[JointDegreeNames.MOTIF_SIZES][i]
                == 0
            )
