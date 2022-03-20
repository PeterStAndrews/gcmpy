import unittest
import numpy as np

from gcmpy.joint_degree.joint_degree_loaders.joint_degree_empirical import JointDegreeEmpirical

NETWORK_SIZE: int = 100000

class JDEmpiricalTest(unittest.TestCase):

    def test_empirical_JDD_single_topology_tuple(self):

        kmin: int = 0
        kmax: int = 25
        n_vertices : int = NETWORK_SIZE
        random_degrees = np.random.randint(kmin,kmax,n_vertices)

        params = {}
        params['jds'] = [(k,) for k in random_degrees]
        params["motif_sizes"] = [2]
                
        DegreeDistObj = JointDegreeEmpirical(params)
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params["motif_sizes"])):
            self.assertTrue(sum([jd[i] for jd in jds]) % params["motif_sizes"][i] == 0)

    def test_empirical_JDD_two_topology_tuple(self):

        kmin: int = 0
        kmax: int = 10
        n_vertices : int = NETWORK_SIZE
        random_degrees = np.random.randint(kmin,kmax,n_vertices)

        params = {}
        params['jds'] = [(k+1,k) for k in random_degrees]
        params["motif_sizes"] = [2,3]
                
        DegreeDistObj = JointDegreeEmpirical(params)
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params["motif_sizes"])):
            self.assertTrue(sum([jd[i] for jd in jds]) % params["motif_sizes"][i] == 0)