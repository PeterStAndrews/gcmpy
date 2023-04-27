import unittest
import networkx as nx

from gcmpy.message_passing.equations.clique_equation import clique_equation
from gcmpy.message_passing.equations.automated_equation import automated_equation
from gcmpy.message_passing.equations.chordless_cycle_equation import chordless_cycle_equation


class EquationTestMixin():

    u: float = 0.651284213
    phi: float = 0.5645231765

    def _initialise_us(self, G: nx.Graph, u: float):
        us = {}
        for n in G.nodes():
            us[n] = u
        nx.set_node_attributes(G, us, "u")

    def _make_chordless_cycle(self, n, u: float) -> nx.Graph:
        G = nx.cycle_graph(n)
        self._initialise_us(G, u)
        return G

    def _make_clique(self, n, u: float) -> nx.Graph:
        G = nx.complete_graph(n)
        self._initialise_us(G, u)
        return G

    def _make_diamond(self, u: float):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
        self._initialise_us(G, u=u)
        return G


class ChordlessCycleEquationTest(EquationTestMixin, unittest.TestCase):
    
    def test_edge(self):
        self.assertAlmostEqual(
            chordless_cycle_equation(3, self.u, self.phi),
            pow(1 - self.phi, 2)
            + 2 * self.phi * pow(1 - self.phi, 2) * self.u
            + (3 * pow(self.phi, 2) * (1 - self.phi) + pow(self.phi, 3))
            * pow(self.u, 2),
            places=7,
        )


class AutomatedEquationTest(EquationTestMixin, unittest.TestCase):

    def test_2_clique(self):
        G = self._make_clique(n=2, u=self.u)
        self.assertAlmostEqual(
            automated_equation(G, self.phi, 0),
            1 - self.phi + self.u * self.phi,
            places=7,
        )

    def test_3_clique(self):
        G = self._make_clique(n=3, u=self.u)
        self.assertAlmostEqual(
            automated_equation(G, self.phi, 0),
            pow(1 - self.phi, 2)
            + 2 * self.phi * pow(1 - self.phi, 2) * self.u
            + (3 * pow(self.phi, 2) * (1 - self.phi) + pow(self.phi, 3))
            * pow(self.u, 2),
            places=7,
        )

    def test_4_cycle(self):
        G = self._make_chordless_cycle(4, self.u)
        self.assertAlmostEqual(
            automated_equation(G, self.phi, 0),
            chordless_cycle_equation(4, self.u, self.phi),
            places=7,
        )

    def test_diamond(self):
        u = self.u
        phi = self.phi

        G = self._make_diamond(u=self.u)

        self.assertAlmostEqual(
            automated_equation(G, self.phi, 0),
            (
                pow(1 - phi, 3)
                + 2 * phi * pow(1 - phi, 3) * u
                + phi * pow(1 - phi, 4) * u
                + 7 * pow(phi, 2) * pow(1 - phi, 3) * pow(u, 2)
                + 2 * pow(phi, 3) * pow(1 - phi, 2) * pow(u, 2)
                + 8 * pow(phi, 3) * pow(1 - phi, 2) * pow(u, 3)
                + 5 * pow(phi, 4) * (1 - phi) * pow(u, 3)
                + pow(phi, 5) * pow(u, 3)
            ),
            places=7,
        )

        self.assertAlmostEqual(
            automated_equation(G, self.phi, 1),
            (
                pow(1 - phi, 2)
                + 2 * phi * pow(1 - phi, 3) * u
                + 5 * pow(phi, 2) * pow(1 - phi, 3) * pow(u, 2)
                + pow(phi, 3) * pow(1 - phi, 2) * pow(u, 2)
                + 8 * pow(phi, 3) * pow(1 - phi, 2) * pow(u, 3)
                + 5 * pow(phi, 4) * (1 - phi) * pow(u, 3)
                + pow(phi, 5) * pow(u, 3)
            ),
            places=7,
        )

    def test_4_clique(self):
        T = self.phi
        u = self.u
        G = self._make_clique(n=4, u=u)
        self.assertAlmostEqual(
            automated_equation(G, self.phi, 1),
            (
                pow(1 - T, 3)
                + 3 * u * T * pow(1 - T, 4)
                + 3
                * pow(u, 2)
                * (pow(T, 3) * pow(1 - T, 3) + 3 * pow(T, 2) * pow(1 - T, 4))
                + pow(u, 3)
                * (
                    pow(T, 6)
                    + 6 * pow(T, 5) * (1 - T)
                    + 15 * pow(T, 4) * pow(1 - T, 2)
                    + 16 * pow(T, 3) * pow(1 - T, 3)
                )
            ),
            places=7,
        )

    def test_n_cycle(self):
        for n in range(3,10):
            G = self._make_chordless_cycle(n, self.u)
            self.assertAlmostEqual(
                automated_equation(G, self.phi, 0),
                chordless_cycle_equation(n, self.u, self.phi),
                places=7
            )

    @unittest.skip("TODO - fix automated_equation() for 5 vertex motifs or more")
    def test_5_clique(self):
        T = self.phi
        u = self.u
        G = self._make_clique(n=5, u=u)

        self.assertAlmostEqual(
            automated_equation(G, T, 0),
            (
                pow(1 - T, 4)
                + 4 * u * T * pow(1 - T, 6)
                + 6
                * pow(u, 2)
                * (pow(T, 3) * pow(1 - T, 6) + 3 * pow(T, 2) * pow(1 - T, 7))
                + 4
                * pow(u, 3)
                * (
                    pow(T, 6) * pow(1 - T, 4)
                    + 6 * pow(T, 5) * pow(1 - T, 5)
                    + 15 * pow(T, 4) * pow(1 - T, 6)
                    + 16 * pow(T, 3) * pow(1 - T, 7)
                )
                + pow(u, 4)
                * (
                    pow(T, 10)
                    + 10 * pow(T, 9) * (1 - T)
                    + 45 * pow(T, 8) * pow(1 - T, 2)
                    + 120 * pow(T, 7) * pow(1 - T, 3)
                    + 205 * pow(T, 6) * pow(1 - T, 4)
                    + 222 * pow(T, 5) * pow(1 - T, 5)
                    + 125 * pow(T, 4) * pow(1 - T, 6)
                )
            ),
            places=7,
        )

    @unittest.skip("TODO - fix automated_equation() for 5 vertex motifs or more")
    def test_6_clique(self):
        T = self.phi
        u = self.u
        G = self._make_clique(n=6, u=u)
        self.assertAlmostEqual(
            automated_equation(G, self.phi, 0),
            pow(1 - T, 5)
            + 5 * u * T * pow(1 - T, 8)
            + 10
            * pow(u, 2)
            * (pow(T, 3) * pow(1 - T, 9) + 3 * pow(T, 2) * pow(1 - T, 10))
            + 10
            * pow(u, 3)
            * (
                pow(T, 6) * pow(1 - T, 8)
                + 6 * pow(T, 5) * pow(1 - T, 9)
                + 15 * pow(T, 4) * pow(1 - T, 10)
                + 16 * pow(T, 3) * pow(1 - T, 11)
            )
            + 5
            * pow(u, 4)
            * (
                pow(T, 10) * pow(1 - T, 5)
                + 10 * pow(T, 9) * pow(1 - T, 6)
                + 45 * pow(T, 8) * pow(1 - T, 7)
                + 120 * pow(T, 7) * pow(1 - T, 8)
                + 205 * pow(T, 6) * pow(1 - T, 9)
                + 222 * pow(T, 5) * pow(1 - T, 10)
                + 125 * pow(T, 4) * pow(1 - T, 11)
            )
            + pow(u, 5)
            * (
                pow(T, 15) + 15 * pow(T, 14) * (1 - T)
                +105 * pow(T, 13) * pow(1 - T, 2)
                + 455 * pow(T, 12) * pow(1 - T, 3)
                + 1365 * pow(T, 11) * pow(1 - T, 4)
                + 2997 * pow(T, 10) * pow(1 - T, 5)
                + 4945 * pow(T, 9) * pow(1 - T, 6)
                + 6165 * pow(T, 8) * pow(1 - T, 7)
                + 5700 * pow(T, 7) * pow(1 - T, 8)
                + 3660 * pow(T, 6) * pow(1 - T, 9)
                + 1296 * pow(T, 5) * pow(1 - T, 10)
            ),
            places=7,
        )


class CliqueEquationTest(EquationTestMixin, unittest.TestCase):

    def test_2_clique(self):
        self.assertAlmostEqual(
            clique_equation(2, self.phi, [self.u]),
            1 - self.phi + self.u * self.phi,
            places=7,
        )

    def test_3_clique(self):
        self.assertAlmostEqual(
            clique_equation(3, self.phi, [self.u, self.u]),
            pow(1 - self.phi, 2)
            + 2 * self.phi * pow(1 - self.phi, 2) * self.u
            + (3 * pow(self.phi, 2) * (1 - self.phi) + pow(self.phi, 3))
            * pow(self.u, 2),
            places=7,
        )

    def test_4_clique(self):
        T = self.phi
        u = self.u
        self.assertAlmostEqual(
            (
                pow(1 - T, 3)
                + 3 * u * T * pow(1 - T, 4)
                + 3
                * pow(u, 2)
                * (pow(T, 3) * pow(1 - T, 3) + 3 * pow(T, 2) * pow(1 - T, 4))
                + pow(u, 3)
                * (
                    pow(T, 6)
                    + 6 * pow(T, 5) * (1 - T)
                    + 15 * pow(T, 4) * pow(1 - T, 2)
                    + 16 * pow(T, 3) * pow(1 - T, 3)
                )
            ),
            clique_equation(4, self.phi, [self.u, self.u, self.u]),
            places=7,
        )

    def test_5_clique(self):
        T = self.phi
        u = self.u
        self.assertAlmostEqual(
            (
                pow(1 - T, 4)
                + 4 * u * T * pow(1 - T, 6)
                + 6
                * pow(u, 2)
                * (pow(T, 3) * pow(1 - T, 6) + 3 * pow(T, 2) * pow(1 - T, 7))
                + 4
                * pow(u, 3)
                * (
                    pow(T, 6) * pow(1 - T, 4)
                    + 6 * pow(T, 5) * pow(1 - T, 5)
                    + 15 * pow(T, 4) * pow(1 - T, 6)
                    + 16 * pow(T, 3) * pow(1 - T, 7)
                )
                + pow(u, 4)
                * (
                    pow(T, 10)
                    + 10 * pow(T, 9) * (1 - T)
                    + 45 * pow(T, 8) * pow(1 - T, 2)
                    + 120 * pow(T, 7) * pow(1 - T, 3)
                    + 205 * pow(T, 6) * pow(1 - T, 4)
                    + 222 * pow(T, 5) * pow(1 - T, 5)
                    + 125 * pow(T, 4) * pow(1 - T, 6)
                )
            ),
            clique_equation(5, self.phi, [self.u] * (5 - 1)),
            places=7,
        )

    def test_6_clique(self):
        T = self.phi
        u = self.u
        self.assertAlmostEqual(
            pow(1 - T, 5)
            + 5 * u * T * pow(1 - T, 8)
            + 10
            * pow(u, 2)
            * (pow(T, 3) * pow(1 - T, 9) + 3 * pow(T, 2) * pow(1 - T, 10))
            + 10
            * pow(u, 3)
            * (
                pow(T, 6) * pow(1 - T, 8)
                + 6 * pow(T, 5) * pow(1 - T, 9)
                + 15 * pow(T, 4) * pow(1 - T, 10)
                + 16 * pow(T, 3) * pow(1 - T, 11)
            )
            + 5
            * pow(u, 4)
            * (
                pow(T, 10) * pow(1 - T, 5)
                + 10 * pow(T, 9) * pow(1 - T, 6)
                + 45 * pow(T, 8) * pow(1 - T, 7)
                + 120 * pow(T, 7) * pow(1 - T, 8)
                + 205 * pow(T, 6) * pow(1 - T, 9)
                + 222 * pow(T, 5) * pow(1 - T, 10)
                + 125 * pow(T, 4) * pow(1 - T, 11)
            )
            + pow(u, 5)
            * (
                pow(T, 15) + 15 * pow(T, 14) * (1 - T)
                +105 * pow(T, 13) * pow(1 - T, 2)
                + 455 * pow(T, 12) * pow(1 - T, 3)
                + 1365 * pow(T, 11) * pow(1 - T, 4)
                + 2997 * pow(T, 10) * pow(1 - T, 5)
                + 4945 * pow(T, 9) * pow(1 - T, 6)
                + 6165 * pow(T, 8) * pow(1 - T, 7)
                + 5700 * pow(T, 7) * pow(1 - T, 8)
                + 3660 * pow(T, 6) * pow(1 - T, 9)
                + 1296 * pow(T, 5) * pow(1 - T, 10)
            ),
            clique_equation(6, T, [u] * (6 - 1)),
            places=7,
        )
