

import unittest
from gcmpy.tools.average_joint_degree_from_jdd import AverageJointDegreeFromJDD

class AverageJointDegreeFromJDDTest(unittest.TestCase):

    def test_average_degrees(self):

        jdd = {(1,2): 0.5, (3,4): 0.5}
        averages = AverageJointDegreeFromJDD.get_average_joint_degrees(jdd)
        expected = [1/2 + 3/2, 2/2 + 4/2]
        self.assertTrue(averages == expected)


        jdd = {(1,2,3): 0.1, (3,4,5): 0.1, (1,1,0): 0.8}
        averages = AverageJointDegreeFromJDD.get_average_joint_degrees(jdd)
        expected = [1 *0.1 + 3 *0.1 + 1*0.8,2 *0.1 + 4 *0.1 + 1*0.8,3 *0.1 + 5 *0.1 + 0*0.8]
        self.assertTrue(averages == expected)
