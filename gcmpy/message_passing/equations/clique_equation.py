import itertools
from gcmpy.message_passing.number_connected_graphs import Q


def clique_equation(tau: int, phi: float, Hs: list) -> float:
    """
    :param tau: size of the clique including focal vertex
    :param phi: bond occupation probability
    :param Hs: H values for all vertices in the clique apart from focal
    """

    def omega(tau, kappa) -> float:
        """
        The number of interface edges for a component of
        kappa vertices in a tau clique.
        :param tau: clique size
        :param kappa: number of neighbours that the focal vertex connects to
        """
        r = tau - kappa - 1
        summation = 0.0
        for v in range(1, r + 1):
            summation += tau - v
        return summation - 0.5 * r * (r - 1)

    summation = 0.0
    # kappa is the number of neighbours that the focal vertex connects to
    for kappa in range(tau):
        factor = []
        for comb in itertools.combinations(Hs, kappa):
            prod = 1
            for H in comb:
                prod *= H

            factor.append(prod)

        for m in range(int(0.5 * kappa * (kappa - 1)) + 1):
            prefactor = (
                Q(kappa + 1, int(0.5 * kappa * (kappa + 1)) - m)
                * pow(phi, int(0.5 * kappa * (kappa + 1)) - m)
                * pow(1 - phi, omega(tau, kappa) + m)
            )

            summation += prefactor * sum(factor)

    return summation
