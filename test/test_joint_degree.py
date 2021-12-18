
# Test the behaviour of the joint degree distribution 
# 
# Copyright (C) 2021 Peter Mann
#
# This file is part of gcmpy, generalised configuration model networks in Python.
#
# gcmpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# gcmpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gcmpy. If not, see <http://www.gnu.org/licenses/gpl.html>.

from gcmpy import *
import unittest

class JDD_manual_Test(unittest.TestCase):

    def test_manual_JDD_single_topology(self):
        
        # valid input data for manual entry
        motif_sizes = [2]
        jdd = {'(1)' : 0.2, '(2)' : 0.5, '(3)' : 0.1, '(5)' : 0.2}

        # create degree distribution
        DegreeDistObj = JDD_manual(jdd, motif_sizes)._jdd

        self.assertEqual(jdd,DegreeDistObj)

    def test_manual_JDD_two_topologies(self):
        
        # valid input data for manual entry
        motif_sizes = [2,3]
        jdd = {'(1,0)' : 0.2, '(2,1)' : 0.5, '(3,0)' : 0.1, '(5,1)' : 0.2}

        # create degree distribution
        DegreeDistObj = JDD_manual(jdd, motif_sizes)._jdd

        self.assertEqual(jdd,DegreeDistObj)




if __name__ == '__main__':
    unittest.main()