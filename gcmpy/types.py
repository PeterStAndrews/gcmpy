
# type hints for gcmpy
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

from typing import NewType, Tuple, Dict, List

_JDD   = NewType('_JDD', Dict[str,float])  
_JDS   = NewType('_JDS', List[List[int]])          
_JOINT_DEGREE = NewType('_JOINT_DEGREE', Tuple[int,...])

_NODE   = NewType('_NODE', int)
_EDGE   = NewType('_EDGE', Tuple[_NODE ,_NODE])    
_NODES  = NewType('_NODES', List[_NODE])
_EDGES  = NewType('_EDGES', List[_EDGE])

_COVER = NewType('_COVER',List[List[int]])

NETWORK_SIZE = 100000
