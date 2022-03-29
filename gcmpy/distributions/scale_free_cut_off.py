import numpy as np


def scale_free_cut_off(alpha: float, kappa: float) -> callable:
    """
    Implements a scale free with exponential degree cutoff function. Limits
    to scale free for large degree cutoff. Undefined for k = 0.
    :param k: int degree
    :param alpha: float power law exponent
    :param kappa: float degree cutoff
    """

    def polylog(s: float, z: float) -> float:
        """
        Implements a polylogarithm function for real arguments, taking two floats
        :param s: base
        :param z: arg
        :returns float: polylogarithm
        """
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

    C = polylog(alpha, np.exp(-1.0 / kappa))  # normalisation constant

    def p(k: int) -> float:
        return (pow((k + 0.0), -alpha) * np.exp(-(k + 0.0) / kappa)) / C

    return p
