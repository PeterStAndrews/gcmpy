import unittest
import networkx as nx

from gcmpy.message_passing.equations.clique_equation import clique_equation
from gcmpy.message_passing.equations.automated_equation import automated_equation


class AutomatedEquationTest(unittest.TestCase):
    u: float = 0.651284213
    phi: float = 0.5645231765

    def _initialise_us(self, G: nx.Graph, u: float = 0.5):
        us = {}
        for n in G.nodes():
            us[n] = u
        nx.set_node_attributes(G, us, "u")

    def _make_edge(self, u: float = 0.5):
        G = nx.Graph()
        G.add_edges_from([(0, 1)])
        self._initialise_us(G, u=u)
        return G

    def _make_triangle(self, u: float = 0.5):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 0)])
        self._initialise_us(G, u=u)
        return G

    def _make_diamond(self, u: float = 0.5):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
        self._initialise_us(G, u=u)
        return G

    def test_2_clique(self):
        G = self._make_edge(u=self.u)
        self.assertAlmostEqual(
            automated_equation(G, self.phi, 0),
            1 - self.phi + self.u * self.phi,
            places=7,
        )

    def test_3_clique(self):
        G = self._make_triangle(u=self.u)
        self.assertAlmostEqual(
            automated_equation(G, self.phi, 0),
            pow(1 - self.phi, 2)
            + 2 * self.phi * pow(1 - self.phi, 2) * self.u
            + (3 * pow(self.phi, 2) * (1 - self.phi) + pow(self.phi, 3))
            * pow(self.u, 2),
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


class CliqueEquationTest(unittest.TestCase):

    u: float = 0.651284213
    phi: float = 0.5645231765

    def _initialise_us(self, G: nx.Graph, u: float = 0.5):
        us = {}
        for n in G.nodes():
            us[n] = u
        nx.set_node_attributes(G, us, "u")

    def test_2_clique(self):
        self.assertAlmostEqual(
            clique_equation(2, self.phi, [self.u]),
            1 - self.phi + self.u * self.phi,
            places=7,
        )

    def test_3_clique(self):
        self.assertAlmostEqual(
            clique_equation(3, self.phi, [self.u,self.u]),
            pow(1 - self.phi, 2)
            + 2 * self.phi * pow(1 - self.phi, 2) * self.u
            + (3 * pow(self.phi, 2) * (1 - self.phi) + pow(self.phi, 3))
            * pow(self.u, 2),
            places=7,
        )

    def test_first_n_cliques(self):
        m = 5
        for n in range(2, m):
            G = nx.complete_graph(n)
            self._initialise_us(G, self.u)
            self.assertAlmostEqual(
                clique_equation(n, self.phi, [self.u]*(n-1)),
                automated_equation(G, self.phi, 0),
                places=7,
            )