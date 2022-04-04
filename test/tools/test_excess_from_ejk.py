import unittest
from gcmpy.tools.joint_excess_from_ejk import JointExcessFromEjk
from gcmpy.tools.joint_excess_joint_degree_matrices import (
    JointExcessJointDegreeMatrices,
)


NETWORK_SIZE: int = 100000


class ExcessFromEjkTest(unittest.TestCase):
    def test_excess_qk_from_ejk(self):

        ejk_tree = {
            (0, 3, 0, 3): 1 / 81,
            (0, 3, 4, 1): 5 / 81,
            (0, 3, 2, 2): 3 / 81,
            (4, 1, 0, 3): 5 / 81,
            (4, 1, 4, 1): 25 / 81,
            (4, 1, 2, 2): 15 / 81,
            (2, 2, 0, 3): 3 / 81,
            (2, 2, 4, 1): 15 / 81,
            (2, 2, 2, 2): 9 / 81,
        }

        ejk_triangle = {
            (3, 1, 3, 1): 16 / 144,
            (3, 1, 1, 2): 24 / 144,
            (3, 1, 5, 0): 8 / 144,
            (1, 2, 3, 1): 24 / 144,
            (1, 2, 1, 2): 36 / 144,
            (1, 2, 5, 0): 12 / 144,
            (5, 0, 3, 1): 8 / 144,
            (5, 0, 1, 2): 12 / 144,
            (5, 0, 5, 0): 4 / 144,
        }

        ejk = JointExcessJointDegreeMatrices()
        ejk._ejks = {"2-clique": ejk_tree, "3-clique": ejk_triangle}
        ejk._excess_degree_keys = {
            "2-clique": [(0, 3), (4, 1), (2, 2)],
            "3-clique": [(3, 1), (1, 2), (5, 0)],
        }

        qks = JointExcessFromEjk.get_excess_joint_distributions(ejk)

        q_tree_test = qks["2-clique"]
        q_triangle_test = qks["3-clique"]

        q_tree_theoretical = {
            (0, 3): 0.1111111111111111,
            (4, 1): 0.5555555555555556,
            (2, 2): 0.3333333333333333,
        }

        q_triangle_theoretical = {
            (3, 1): 0.33333333333333337,
            (1, 2): 0.49999999999999994,
            (5, 0): 0.16666666666666669,
        }

        self.assertTrue(q_tree_test == q_tree_theoretical)
        self.assertTrue(q_triangle_test, q_triangle_theoretical)
