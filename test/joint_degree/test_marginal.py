import unittest

from gcmpy.joint_degree.joint_degree_loaders.joint_degree_marginal import (
    JointDegreeMarginal,
)
from gcmpy.names.joint_degree_names import JointDegreeNames
from gcmpy.distributions.poisson import poisson
from gcmpy.distributions.power_law import power_law
from gcmpy.distributions.scale_free_cut_off import scale_free_cut_off

NETWORK_SIZE: int = 100000


class JDMarginalTest(unittest.TestCase):
    def test_marginal_JDD_single_topology(self):

        n_vertices: int = NETWORK_SIZE
        mean_degree: float = 2.5  # mean poisson degree
        kmin: int = 0  # smallest degree
        kmax: int = 10  # largest degree

        params = {}
        params[JointDegreeNames.MOTIF_SIZES] = [2]
        params[JointDegreeNames.ARR_FP] = [poisson(mean_degree)]
        params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND] = [(kmin, kmax)]

        DegreeDistObj = JointDegreeMarginal(params)
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params[JointDegreeNames.MOTIF_SIZES])):
            self.assertTrue(
                sum([jd[i] for jd in jds]) % params[JointDegreeNames.MOTIF_SIZES][i]
                == 0
            )

    def test_marginal_JDD_two_topologies(self):

        n_vertices: int = NETWORK_SIZE
        mean_degree: float = 2.5  # mean poisson degree
        kmin: int = 0  # smallest degree
        kmax: int = 10  # largest degree

        params = {}
        params[JointDegreeNames.MOTIF_SIZES] = [2, 3]
        params[JointDegreeNames.ARR_FP] = [poisson(mean_degree), poisson(mean_degree)]

        params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND] = [(kmin, kmax)] * len(
            params[JointDegreeNames.MOTIF_SIZES]
        )

        DegreeDistObj = JointDegreeMarginal(params)
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params[JointDegreeNames.MOTIF_SIZES])):
            self.assertTrue(
                sum([jd[i] for jd in jds]) % params[JointDegreeNames.MOTIF_SIZES][i]
                == 0
            )

    def test_marginal_JDD_two_topologies_sampling(self):

        n_vertices: int = NETWORK_SIZE
        mean_degree: float = 2.5  # mean poisson degree
        kmin: int = 0  # smallest degree
        kmax: int = 10  # largest degree

        params = {}
        params[JointDegreeNames.MOTIF_SIZES] = [2, 3]
        params[JointDegreeNames.ARR_FP] = [poisson(mean_degree), poisson(mean_degree)]

        params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND] = [(kmin, kmax)] * len(
            params[JointDegreeNames.MOTIF_SIZES]
        )
        params[JointDegreeNames.USE_SAMPLING] = True
        params[JointDegreeNames.N_SAMPLES] = 100000

        DegreeDistObj = JointDegreeMarginal(params)
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params[JointDegreeNames.MOTIF_SIZES])):
            self.assertTrue(
                sum([jd[i] for jd in jds]) % params[JointDegreeNames.MOTIF_SIZES][i]
                == 0
            )

    def test_marginal_JDD_power_law(self):

        n_vertices: int = NETWORK_SIZE
        powerlaw_exponent = 2.5  # powerlaw_exponent
        kmin: int = 1  # smallest degree
        kmax: int = 10  # largest degree

        params = {}
        params[JointDegreeNames.MOTIF_SIZES] = [2]
        params[JointDegreeNames.ARR_FP] = [power_law(powerlaw_exponent)]
        params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND] = [(kmin, kmax)]

        DegreeDistObj = JointDegreeMarginal(params)
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params[JointDegreeNames.MOTIF_SIZES])):
            self.assertTrue(
                sum([jd[i] for jd in jds]) % params[JointDegreeNames.MOTIF_SIZES][i]
                == 0
            )

    def test_marginal_JDD_scale_free(self):

        n_vertices: int = NETWORK_SIZE
        powerlaw_exponent = 2.5  # powerlaw_exponent
        degree_cut_off = 25  # degree cut off
        kmin: int = 1  # smallest degree
        kmax: int = 10  # largest degree

        params = {}
        params[JointDegreeNames.MOTIF_SIZES] = [2]
        params[JointDegreeNames.LOW_HIGH_DEGREE_BOUND] = [(kmin, kmax)]
        params[JointDegreeNames.ARR_FP] = [
            scale_free_cut_off(powerlaw_exponent, degree_cut_off)
        ]

        DegreeDistObj = JointDegreeMarginal(params)
        jds = DegreeDistObj.sample_jds_from_jdd(n_vertices)

        self.assertTrue(len(jds) == n_vertices)
        for i in range(len(params[JointDegreeNames.MOTIF_SIZES])):
            self.assertTrue(
                sum([jd[i] for jd in jds]) % params[JointDegreeNames.MOTIF_SIZES][i]
                == 0
            )
