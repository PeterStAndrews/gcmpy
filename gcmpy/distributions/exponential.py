import numpy as np


def exponential(a: float) -> callable:
    """
    Implemnts an exponential distribution.
    :param a: distribution parameter, a>0.
    :returns p: callable
    """

    def p(k: int) -> float:
        return (1 - np.exp(-a)) * np.exp(-a * k)

    return p
