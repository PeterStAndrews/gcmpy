
from gcmpy import *
import unittest

def EECC_Test_util_network()->EECC:

    G = EECC()

    G.add_edge((1, 2))
    G.add_edge((1, 14))
    G.add_edge((2, 4))
    G.add_edge((2, 13))
    G.add_edge((2, 14))
    G.add_edge((3, 4))
    G.add_edge((3, 5))
    G.add_edge((4, 5))
    G.add_edge((4, 13))
    G.add_edge((4, 14))
    G.add_edge((6, 7))
    G.add_edge((6, 13))
    G.add_edge((7, 8))
    G.add_edge((7, 13))
    G.add_edge((8, 9))
    G.add_edge((8, 13))
    G.add_edge((9, 10))
    G.add_edge((9, 11))
    G.add_edge((9, 13))
    G.add_edge((10, 11))
    G.add_edge((11, 12))
    G.add_edge((12, 13))
    G.add_edge((13, 14))

    return G
        
class EECC_Test(unittest.TestCase):

    def test_EECC_m0_2(self):
        '''Test adapted from (2) for python.

         ----- References -----

        1) Giulio Burgio, Alex Arenas, Sergio Gómez and Joan T. Matamalas: 
            Network clique cover approximation to analyze complex contagions through group interactions, 
            Comms. Phys. (2021) in press (arXiv:2101.03618)

        2) https://github.com/giubuig/DisjointCliqueCover.jl'''

        
        G = EECC_Test_util_network()

        # Sets of maximal cliques in G to test against
        C_2 = [[1, 2], [1, 14], [2, 4], [2, 13], [2, 14], [3, 4], [3, 5], [4, 5], [4, 13],
               [4, 14], [6, 7], [6, 13], [7, 8], [7, 13], [8, 9], [8, 13], [9, 10], [9, 11],
                [9, 13], [10, 11], [11, 12], [12, 13], [13, 14]]

        # Cliques extraction

        G.set_max_clique_size(2)
        C_2_test = G.limited_maximal_cliques()
        self.assertEqual(C_2,C_2_test)

        # scores of the maximal cliques in G
        r_2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        EC_2 : List[int] = []
        ord_2 : List[int] = [0] * len(C_2)
        r_2_test : List[float] = [0.0] * len(C_2)
        indexes_2 : List[int] = []
        
        G.set_max_clique_size(2)
        G.compute_scores(C_2, EC_2, ord_2, r_2_test, indexes_2)

        for i, j in zip(r_2,r_2_test):
            self.assertAlmostEqual(i, j, 1e-5)

        EECC_2 = [[1, 2], [1, 14], [2, 4], [2, 13], [2, 14], [3, 4], [3, 5], [4, 5], [4, 13],
                [4, 14], [6, 7], [6, 13], [7, 8], [7, 13], [8, 9], [8, 13], [9, 10], [9, 11],
                [9, 13], [10, 11], [11, 12], [12, 13], [13, 14]]

        # EECC calculation
        G.set_max_clique_size(2)
        EECC_2_test = G.get_EECC()

        self.assertEqual(EECC_2,EECC_2_test)
       
        

    def test_EECC_m0_3(self):

        G = EECC_Test_util_network()

        C_3 = [[1, 2, 14], [2, 4, 13], [2, 4, 14], [2, 13, 14], [3, 4, 5], [4, 13, 14],
                [6, 7, 13], [7, 8, 13], [8, 9, 13], [9, 10, 11], [11, 12], [12, 13]]

        G.set_max_clique_size(3)
        C_3_test = G.limited_maximal_cliques()
        self.assertEqual(C_3,C_3_test)

        r_3 = [1/3, 1.0, 1.0, 1.0, 0.0, 1.0, 1/3, 2/3, 1/3, 0.0, 0.0, 0.0]

        EC_3 : List[int] = []
        ord_3 : List[int] = [0] * len(C_3)
        r_3_test : List[float] = [0.0] * len(C_3)
        indexes_3 : List[int] = []
       
        G.compute_scores(C_3, EC_3, ord_3, r_3_test, indexes_3)
        
        for i, j in zip(r_3,r_3_test):
            self.assertAlmostEqual(i, j, 1e-5)

        EECC_3 = [[1, 2, 14], [3, 4, 5], [4, 13, 14], [6, 7, 13], [8, 9, 13], [9, 10, 11],
                 [2, 4], [2, 13], [7, 8], [11, 12], [12, 13]]

       
        EECC_3_test = G.get_EECC()

        self.assertEqual(EECC_3,EECC_3_test)

    def test_EECC_m0_4(self):

        G = EECC_Test_util_network()

        C_4 = [[2, 4, 13, 14], [1, 2, 14], [3, 4, 5], [6, 7, 13], [7, 8, 13], [8, 9, 13],
                [9, 10, 11], [11, 12], [12, 13]]

        G.set_max_clique_size(4)
        C_4_test = G.limited_maximal_cliques()
        self.assertEqual(C_4,C_4_test)

        r_4 = [1/6, 1/3, 0.0, 1/3, 2/3, 1/3, 0.0, 0.0, 0.0]

        EC_4 : List[int] = []
        ord_4 : List[int] = [0] * len(C_4)
        r_4_test : List[float] = [0.0] * len(C_4)
        indexes_4 : List[int] = []
        G.set_max_clique_size(4)
        G.compute_scores(C_4, EC_4, ord_4, r_4_test, indexes_4)
        
        for i, j in zip(r_4,r_4_test):
            self.assertAlmostEqual(i, j, 1e-5)

        EECC_4 = [[2, 4, 13, 14], [3, 4, 5], [6, 7, 13], [8, 9, 13], [9, 10, 11],
                    [1, 2], [1, 14], [7, 8], [11, 12], [12, 13]]

        G.set_max_clique_size(4)
        EECC_4_test = G.get_EECC()

        self.assertEqual(EECC_4,EECC_4_test)