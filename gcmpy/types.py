
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

from typing import Tuple, TypeVar, Dict, List
from gcm_algorithm import output_data

_JDD   = TypeVar('_JDD', bound=Dict[str,float])  
_JDS   = TypeVar('_JDS', bound=List[List[int]])          
_JOINT_DEGREE = TypeVar('_JOINT_DEGREE', bound=Tuple[int,...])

_EDGE   = TypeVar('_EDGE', bound=Tuple[int,int])    
_NODES  = TypeVar('_NODES', bound=List[int])

_RESULTS   = TypeVar('_RESULTS', bound=List[output_data])  
