# Utilities for gcmpy
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

import io
import pickle
from typing import List

from .types import _EDGES

class edge_list(object):
    '''Network represented as an edge list. The GCM class uses this 
    structure to generate networks which can then be converted to 
    other network libraries.'''
    
    def __init__(self):
        self._edge_list : _EDGES = []
        
    def add_edges_from(self, edges : _EDGES)->None:
        '''Adds edges from list of tuples (int,int) to the edge list.
        :param edges: list of tuples of ints.'''
        for e in edges:
            self._edge_list.append(e)

class output_data(object):
    '''An object to store output data from the process.
    :param i: integer for experiment index'''

    def __init__(self, i : int):
        self._experiment : int = i                # experiment index
        self._name : str = ''                     # tags for network
        self._network  : edge_list = None         # network
    
class results(object):
    '''A collection of output_data objects that can
    be serialised and converted to other graph formats.'''

    def __init__(self):
        self.res : List[output_data]

    def add_result(self, r : output_data)->None:
        self.res.append(r)

    def serialise_results_to_file(self, filename : str)->None:
        '''Dump the results structure to a binary file.
        :param filename: name of file to create.'''
        
        binary_file = open(filename, mode='wb')
        pickle.dump(self.res, binary_file)
        binary_file.close()

    def read_results_from_binary_file(self, filename : str)->None:
        '''read the results structure from a binary file.
        :param filename: name of the file to read.'''

        f : io.BufferedReader = open(filename,mode='rb')
        bin_data : bytes = f.read()
        sio : io.StringIO = io.StringIO(bin_data)
        self.res : results = pickle.load(sio)

    def to_networkx(self)->List[any]:
        pass

    def to_iGraph(self)->List[any]:
        pass