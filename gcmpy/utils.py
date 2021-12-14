

from typing import Callable
import numpy as np

def poisson(kmean : float)->Callable:
    '''Implements a poisson distribution.'''
    def p(k : int)->float:
        return np.exp() * pow(kmean,k) / np.factorial(k)
    return p

def exponential(a : float)->Callable:
    '''Implemnts an exponential distribution.'''
    def p(k : int)->float:
        return (1-a)*pow(a,k)
    return p

def power_law(alpha : float)->Callable:
    '''Implements a power law distribution with exponent alpha. Undefined for k = 0.'''

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
    Param k: degree
    Param alpha: exponent
    Param kappa: degree cutoff'''

    def polylog(s : float, z : float) -> float:
        '''Implements a polylogarithm function for real arguments, taking two floats
        param s: base
        param z: arg
        returns polylogarithm float.'''
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

