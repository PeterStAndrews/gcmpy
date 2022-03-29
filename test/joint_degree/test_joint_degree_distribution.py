
import unittest
import numpy as np

from gcmpy.joint_degree.joint_degree import JointDegree
from gcmpy.joint_degree.joint_degree_types import JointDegreeType
from gcmpy.joint_degree.joint_degree_distribution import JointDegreeDistribution
from gcmpy.names.joint_degree_names import JointDegreeNames
from gcmpy.distributions.poisson import poisson
from gcmpy.distributions.power_law import power_law

NETWORK_SIZE: int = 100000

class JointDegreeDistributionTest(unittest.TestCase):

    def test_manual(self):

        params = {}
        params[JointDegreeNames.JOINT_DEGREE_TYPE] = JointDegreeType.MANUAL
        params[JointDegreeNames.JDD] = {(1,) : 0.2, (2,) : 0.5, (3,) : 0.1, (5,) : 0.2}
        params[JointDegreeNames.MOTIF_SIZES] = [2]

        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)
        
        n_vertices : int = NETWORK_SIZE
        jds = DegreeDist.sample_jds_from_jdd(n_vertices)
           
        
    def test_empirical(self):

        kmin: int = 0
        kmax: int = 25
        n_vertices : int = NETWORK_SIZE
        random_degrees = np.random.randint(kmin,kmax,n_vertices)

        params = {}
        params[JointDegreeNames.JOINT_DEGREE_TYPE] = JointDegreeType.EMPIRICAL
        params[JointDegreeNames.JDS] = [(k,) for k in random_degrees]
        params[JointDegreeNames.MOTIF_SIZES] = [2]
                
        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)

        jds = DegreeDist.sample_jds_from_jdd(n_vertices)

    def test_marginal(self):

        n_vertices: int = NETWORK_SIZE
        mean_degree: float = 2.5                    # mean poisson degree
        kmin: int = 0                               # smallest degree
        kmax: int = 10                              # largest degree

        params = {}
        params[JointDegreeNames.JOINT_DEGREE_TYPE] = JointDegreeType.MARGINAL
        params[JointDegreeNames.MOTIF_SIZES] = [2]
        params[JointDegreeNames.ARR_FP] = [poisson(mean_degree)] 
        params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND] = [(kmin,kmax)]

        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)

        jds = DegreeDist.sample_jds_from_jdd(n_vertices)

    def test_split_degree(self):
     
        n_vertices: int = NETWORK_SIZE
        kmax: int = 1000                                  
        kmin: int = 1  
        powerlaw_exponent: float = 2.5   

        params = {}
        params[JointDegreeNames.JOINT_DEGREE_TYPE] = JointDegreeType.SPLIT_DEGREE
        params[JointDegreeNames.MOTIF_SIZES] = [2,3]                    
        params[JointDegreeNames.PROBS] = [0.8,0.2]                       # prob edge is 2- or 3-clique
        params[JointDegreeNames.FP] = power_law(powerlaw_exponent)       # overall degree distribution     
        params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND] = (kmin,kmax)

        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)

        jds = DegreeDist.sample_jds_from_jdd(n_vertices)

    def test_delta(self):

        n_vertices: int = NETWORK_SIZE
        kmax: int = 1000                                  
        kmin: int = 1  
        powerlaw_exponent: float = 2.5   

        params = {}
        params[JointDegreeNames.JOINT_DEGREE_TYPE] = JointDegreeType.DELTA
        params[JointDegreeNames.MOTIF_SIZES] = [2,3]                    
        params[JointDegreeNames.PROBS] = [0.8,0.2]                       # prob edge is 2- or 3-clique
        params[JointDegreeNames.FP] = power_law(powerlaw_exponent)       # overall degree distribution     
        params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND] = (kmin,kmax)
        params[JointDegreeNames.TARGET_K] = 3
        
        DegreeDist: JointDegree = JointDegreeDistribution.load_joint_degree(params)

        jds = DegreeDist.sample_jds_from_jdd(n_vertices)