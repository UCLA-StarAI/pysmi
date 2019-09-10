import unittest

import networkx as nx
from pysmt.shortcuts import Symbol, Real, And, LE, Or
from pysmt.typing import REAL

from pysmi.graph2formula import check_tree, graph_to_tree
from pysmi.smtree import TNode, TEdge
from pysmi.smi import SMI


class TestSMI(unittest.TestCase):
    def test_pe_edge(self):
        sym_x = Symbol("x", REAL)  # child
        sym_y = Symbol("y", REAL)  # parent

        x_domain = And(sym_x > Real(-1), sym_x < Real(6))
        y_domain = And(sym_y > Real(-2), sym_y < Real(4))

        x = TNode(symbol=sym_x, label=1, domains=x_domain)
        y = TNode(symbol=sym_y, label=0, domains=y_domain)

        formula = LE(sym_x, sym_y + Real(1))
        edge = TEdge(parent=y, child=x, formula=formula)
        y.add_edge(edge)

        x_constraints = [[0, 1, 2], [2.5, 3.5, 4]]
        interval_and_degree = SMI.pe_edge(edge, x_constraints)
        self.assertListEqual(interval_and_degree[0], [-1.0, 0.0, 3])
        self.assertListEqual(interval_and_degree[1], [0, 1.5, 0])
        self.assertListEqual(interval_and_degree[2], [1.5, 2.5, 5])
        self.assertListEqual(interval_and_degree[3], [2.5, 4, 0])

    def test_pe_node(self):
        sym_x0 = Symbol("x0", REAL)
        sym_x1 = Symbol("x1", REAL)
        sym_x2 = Symbol("x2", REAL)

        formula1 = Or(sym_x0 + Real(1) < sym_x1, sym_x0 - Real(1) > sym_x1)
        formula2 = Or(sym_x1 + Real(1) < sym_x2, sym_x1 - Real(1) > sym_x2)

        x0 = TNode(sym_x0, 0, n_bbox=4, xrange=5)
        x1 = TNode(sym_x1, 1, n_bbox=4, xrange=5)
        x2 = TNode(sym_x2, 2)

        edge1 = TEdge(x0, x1, formula1)
        x0.add_edge(edge1)
        edge2 = TEdge(x1, x2, formula2)
        x1.add_edge(edge2)

        check_tree(x0)

        interval_and_degree = SMI.pe_node(x0)

        # self.assertListEqual(interval_and_degree[0], [-1.0, 0.0, 2])
        # self.assertListEqual(interval_and_degree[1], [0, 1, 2])

    def test_compute_mi_path(self):
        gold = [2.0, 1.0, 0.6666666666, 0.4166666666]
        graphs = [nx.generators.path_graph(i) for i in range(1, 1 + len(gold))]
        smi = SMI()
        res = [smi.compute_mi(graph_to_tree(g, n_bbox=1)) for g in graphs]
        for a, b in zip(res, gold):
            self.assertAlmostEqual(a, b)

    def test_compute_mi_fat(self):
        gold = [2.0, 1.0, 0.6666666666666667, 0.5000000000000011,
                0.2999999999999875, 0.21111111111116665]
        branches = 3
        graphs = [nx.full_rary_tree(branches, i)
                  for i in range(1, 1 + len(gold))]
        smi = SMI()
        res = [smi.compute_mi(graph_to_tree(g, n_bbox=1)) for g in graphs]
        for a, b in zip(res, gold):
            self.assertAlmostEqual(a, b)

    def test_compute_mi_star(self):
        gold = [1.0, 0.6666666666666667, 0.5000000000000011,
                0.4000000000000046, 0.33333333333332393, 0.28571428571429613]
        graphs = [nx.star_graph(i) for i in range(1, 1 + len(gold))]
        # res = [SMI.compute_mi(graph_to_tree(g)) for g in graphs]
        smi = SMI()
        res = [smi.compute_mi(graph_to_tree(g, n_bbox=1)) for g in graphs]
        print("res:{}".format(res))
        for a, b in zip(res, gold):
            self.assertAlmostEqual(a, b)
