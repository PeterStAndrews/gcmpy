

import unittest
from gcmpy.tools.joint_excess_from_joint_degree import JointExcessDDfromJDD

class JointExcessDDfromJDDTest(unittest.TestCase):

    def test_joint_excess_from_jdd(self):

        jdd = {(5,1) : 1/3, (3,2) : 1/3, (1,3) : 1/3}

        q_tree_st_expected = {(0, 3): 0.1111111111111111,
                              (4, 1): 0.5555555555555556,
                              (2, 2): 0.3333333333333333}

        q_triangle_st_expected = {(3, 1): 0.33333333333333337,
                                  (1, 2): 0.49999999999999994,
                                  (5, 0): 0.16666666666666669}

        test_data:list[dict] = JointExcessDDfromJDD.get_joint_excess_distributions(jdd)

        for key in q_tree_st_expected:
            self.assertAlmostEqual(q_tree_st_expected[key], test_data[0][key])

        for key in q_triangle_st_expected:
            self.assertAlmostEqual(q_triangle_st_expected[key], test_data[1][key])

        test_data_dict:dict[dict] = JointExcessDDfromJDD.convert_list_qks_to_dict(
            test_data,['2-clique', '3-clique']
        )

        self.assertTrue(test_data_dict['2-clique'] == test_data[0])
        self.assertTrue(test_data_dict['3-clique'] == test_data[1])