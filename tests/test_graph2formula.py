import unittest

from pysmt.shortcuts import Real, Symbol, Or
from pysmt.typing import REAL

from pysmi.graph2formula import tree_to_formula
from pysmi.smtree import TNode, TEdge


class TestGraph2Formula(unittest.TestCase):
    def test_tree_to_formula(self):
        sym_x0 = Symbol("x0", REAL)
        sym_x1 = Symbol("x1", REAL)
        sym_x2 = Symbol("x2", REAL)

        formula1 = Or(sym_x0 + Real(1) < sym_x1, sym_x0 - Real(1) > sym_x1)
        formula2 = Or(sym_x1 + Real(1) < sym_x2, sym_x1 - Real(1) > sym_x2)

        x0 = TNode(sym_x0, 0)
        x1 = TNode(sym_x1, 1)
        x2 = TNode(sym_x2, 2)

        edge1 = TEdge(x0, x1, formula1)
        x0.add_edge(edge1)
        edge2 = TEdge(x1, x2, formula2)
        x1.add_edge(edge2)

        formula = tree_to_formula(x0)
        # TODO
