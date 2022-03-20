
import unittest
import numpy as np

from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_types import JointDegreeType
from gcmpy.joint_degree.joint_degree_distribution import JointDegreeDistribution

from gcmpy.distributions.poisson import poisson
from gcmpy.distributions.power_law import power_law

NETWORK_SIZE: int = 100000

class JointDegreeDistributionTest(unittest.TestCase):

    def test_manual(self):

        params = {}
        params["joint_degree_type"] = JointDegreeType.MANUAL
        params["jdd"] = {(1,) : 0.2, (2,) : 0.5, (3,) : 0.1, (5,) : 0.2}
        params["motif_sizes"] = [2]

        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)
        
        n_vertices : int = NETWORK_SIZE
        jds = DegreeDist.sample_jds_from_jdd(n_vertices)
           
        
    def test_empirical(self):

        kmin: int = 0
        kmax: int = 25
        n_vertices : int = NETWORK_SIZE
        random_degrees = np.random.randint(kmin,kmax,n_vertices)

        params = {}
        params["joint_degree_type"] = JointDegreeType.EMPIRICAL
        params['jds'] = [(k,) for k in random_degrees]
        params["motif_sizes"] = [2]
                
        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)

        jds = DegreeDist.sample_jds_from_jdd(n_vertices)

    def test_marginal(self):

        n_vertices: int = NETWORK_SIZE
        mean_degree: float = 2.5                    # mean poisson degree
        kmin: int = 0                               # smallest degree
        kmax: int = 10                              # largest degree

        params = {}
        params["joint_degree_type"] = JointDegreeType.MARGINAL
        params["motif_sizes"] = [2]
        params["arr_fp"] = [poisson(mean_degree)] 
        params["low_high_degree_bounds"] = [(kmin,kmax)]

        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)

        jds = DegreeDist.sample_jds_from_jdd(n_vertices)

    def test_split_degree(self):
     
        n_vertices: int = NETWORK_SIZE
        kmax: int = 1000                                  
        kmin: int = 1  
        powerlaw_exponent: float = 2.5   

        params = {}
        params["joint_degree_type"] = JointDegreeType.SPLIT_DEGREE
        params["motif_sizes"] = [2,3]                    
        params["probs"] = [0.8,0.2]                       # prob edge is 2- or 3-clique
        params["fp"] = power_law(powerlaw_exponent)       # overall degree distribution     
        params["low_high_degree_bound"] = (kmin,kmax)

        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)

        jds = DegreeDist.sample_jds_from_jdd(n_vertices)

    def test_delta(self):

        n_vertices: int = NETWORK_SIZE
        kmax: int = 1000                                  
        kmin: int = 1  
        powerlaw_exponent: float = 2.5   

        params = {}
        params["joint_degree_type"] = JointDegreeType.DELTA
        params["motif_sizes"] = [2,3]                    
        params["probs"] = [0.8,0.2]                       # prob edge is 2- or 3-clique
        params["fp"] = power_law(powerlaw_exponent)       # overall degree distribution     
        params["low_high_degree_bound"] = (kmin,kmax)
        params["target_k"] = 3
        
        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)

        jds = DegreeDist.sample_jds_from_jdd(n_vertices)