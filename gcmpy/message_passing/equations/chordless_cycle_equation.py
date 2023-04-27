def chordless_cycle_equation(n, u, phi):
    """
    P. Mann (et al). Phys. Rev. E 103, 012313 2021
    """
    summation = sum(
        [(i + 1) * pow(phi * u, i) * pow(1 - phi, 2) for i in range(1, n - 1)]
    )
    return (
        pow(1 - phi, 2)
        + summation
        + n * pow(u * phi, n - 1) * (1 - phi)
        + phi * pow(phi * u, n - 1)
    )
