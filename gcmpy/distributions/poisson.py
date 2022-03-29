import numpy as np


def poisson(kmean: float) -> callable:
    """
    Implements a poisson distribution
    :param kmean: mean of poisson distribution
    :returns p: Callable
    """

    def p(k: int) -> float:
        return np.exp(-kmean) * pow(kmean, k) / np.math.factorial(k)

    return p
