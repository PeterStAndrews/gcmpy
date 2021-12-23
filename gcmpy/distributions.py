# distribution functions for gcmpy
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

from typing import Callable
import numpy as np

def poisson(kmean : float)->Callable:
    '''Implements a poisson distribution.
    
    :param kmean: mean of poisson distribution
    :returns p: Callable'''
    def p(k : int)->float:
        return np.exp(-kmean) * pow(kmean,k) / np.math.factorial(k)
    return p

def exponential(a : float)->Callable:
    '''Implemnts an exponential distribution.
    :param a: distribution parameter, a>0.
    :returns p: callable'''
    def p(k : int)->float:
        return (1-np.exp(-a))*np.exp(-a*k)
    return p

def power_law(alpha : float)->Callable:
    '''Implements a power law distribution with exponent alpha. Undefined for k = 0.
    
    :param alpha: power law exponent
    :returns p: callable'''

    def zeta(s : float)->float:
        tol = +1e-06
        l = 0.0
        k = 1
        while 1:
            term = 1.0/ k**s
            l += term
            if abs(term) < tol:
                break
            k += 1
        return l

    C = zeta(alpha)

    def p(k:int)->float:
        return pow(k,-alpha) / C
    return p

def scale_free_cut_off(alpha : float, kappa :float) -> Callable:
    '''Implements a scale free with exponential degree cutoff function. Limits
    to scale free for large degree cutoff. Undefined for k = 0.
    
    :param k: int degree
    :param alpha: float power law exponent
    :param kappa: float degree cutoff'''

    def polylog(s : float, z : float) -> float:
        '''Implements a polylogarithm function for real arguments, taking two floats
        
        :param s: base
        :param z: arg
        :returns polylogarithm float:'''
        tol = +1e-06
        l = 0.0
        k = 1
        zk = z
        while 1:
            term = zk / k**s
            l += term
            if abs(term) < tol:
                break
            zk *= z
            k += 1
        return l

    C = polylog(alpha, np.exp(-1.0 / kappa)) # normalisation constant

    def p(k : int) -> float:
        return (pow((k + 0.0), -alpha) * np.exp(-(k + 0.0) / kappa)) / C
    return p

def log_normal():
    pass

def weibull():
    pass