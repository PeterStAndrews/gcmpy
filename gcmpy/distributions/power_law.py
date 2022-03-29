def power_law(alpha: float) -> callable:
    """
    Implements a power law distribution with exponent alpha. Undefined for k = 0.
    :param alpha: power law exponent
    :returns p: callable
    """

    def zeta(s: float) -> float:
        tol = +1e-06
        l = 0.0
        k = 1
        while 1:
            term = 1.0 / k**s
            l += term
            if abs(term) < tol:
                break
            k += 1
        return l

    C = zeta(alpha)

    def p(k: int) -> float:
        return pow(k, -alpha) / C

    return p
