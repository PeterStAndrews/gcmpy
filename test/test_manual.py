
import unittest

from gcmpy.joint_degree.joint_degree_loaders.joint_degree_manual import JointDegreeManual

NETWORK_SIZE: int = 100000

class JDManualTest(unittest.TestCase):

    def test_manual_JDD_single_topology(self):

        # valid input data for manual entry
        params = {}
        params["jdd"] = {(1,) : 0.2, (2,) : 0.5, (3,) : 0.1, (5,) : 0.2}
        params["motif_sizes"] = [2]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices : int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params["motif_sizes"])):
            self.assertTrue(
                sum([jd[0] for jd in jds]) % params["motif_sizes"][i] == 0
            )

    def test_manual_JDD_two_topologies(self):
        
        params = {}
        params["jdd"] = {(1,0) : 0.2, (2,1) : 0.5, (3,0) : 0.1, (5,1) : 0.2}
        params["motif_sizes"] = [2,3]
        
        DegreeDistObj = JointDegreeManual(params)
        n_vertices : int = NETWORK_SIZE 
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params["motif_sizes"])):
            self.assertTrue(
                sum([jd[i] for jd in jds]) % params["motif_sizes"][i] == 0
                )