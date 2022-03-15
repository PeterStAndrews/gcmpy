
import unittest

from gcmpy.joint_degree.joint_degree_loaders.joint_degree_split_degree import JointDegreeSplitDegree
from gcmpy.distributions.power_law import power_law

NETWORK_SIZE: int = 100000

class JDSplitDegreeTest(unittest.TestCase):

    def test_split_k_two_topologies(self):

        kmax: int = 1000                                  
        kmin: int = 1  
        powerlaw_exponent: float = 2.5   

        params = {}
        params["motif_sizes"] = [2,3]                    
        params["probs"] = [0.8,0.2]                       # prob edge is 2- or 3-clique
        params["fp"] = power_law(powerlaw_exponent)       # overall degree distribution     
        params["low_high_degree_bound"] = (kmin,kmax)

        DegreeDistObj = JointDegreeSplitDegree(params)

        n_vertices: int = NETWORK_SIZE
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params["motif_sizes"])):
            self.assertTrue(sum([jd[i] for jd in jds]) % params["motif_sizes"][i] == 0)


