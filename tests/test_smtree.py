import unittest

from pysmt.shortcuts import Symbol, Real, Or, LT, LE, GE, And
from pysmt.typing import REAL

from pysmi.smtree import TNode, TEdge
from pysmi.utils import get_constants, get_x_coefficient, \
    get_real_variables, categorize_bounds


class TestSmtree(unittest.TestCase):
    def setUp(self) -> None:
        parent_symbol = Symbol("x", REAL)
        child_symbol = Symbol("y", REAL)
        # y < 2x+1 or y > x - 1 or x > 2 or y > 0.1
        atom0 = LE(child_symbol, Real(2) * parent_symbol + Real(1))
        atom1 = LE(parent_symbol - 1, child_symbol)
        atom2 = LE(Real(2), parent_symbol)
        atom3 = LE(Real(0.1), child_symbol)
        atom4 = LE(child_symbol, Real(0.4))
        self.parent_symbol = parent_symbol
        self.child_symbol = child_symbol
        self.atoms = [atom0, atom1, atom2, atom3, atom4]
        self.formula = Or(self.atoms[:4])

    def test_categorize_bounds(self):
        clause = self.formula
        formula = And(clause, clause, clause)

        upper, lower, parent_atoms = categorize_bounds(
            formula,
            self.child_symbol)
        # upper
        const = get_constants(upper[0])
        coef = get_x_coefficient(upper[0])
        self.assertAlmostEqual(float(const.constant_value()), 1)
        self.assertAlmostEqual(float(coef.constant_value()), 2)
        # lower
        l_idx = 0 if len(get_real_variables(lower[0])) != 0 else 1
        r_idx = abs(1 - l_idx)

        const = get_constants(lower[l_idx])
        coef = get_x_coefficient(lower[l_idx])
        self.assertAlmostEqual(float(const.constant_value()), -1)
        self.assertAlmostEqual(float(coef.constant_value()), 1)
        const = get_constants(lower[r_idx])
        self.assertAlmostEqual(float(const.constant_value()), 0.1)
        # parent atoms
        self.assertListEqual(parent_atoms, [self.atoms[2]])

        # another test
        formula = And(self.atoms[3], self.atoms[4])
        upper, lower, parent_atoms = categorize_bounds(
            formula,
            self.child_symbol)
        # upper
        const = get_constants(upper[0])
        self.assertAlmostEqual(float(const.constant_value()), 0.4)
        # lower
        const = get_constants(lower[0])
        self.assertAlmostEqual(float(const.constant_value()), 0.1)

    def test_collect_bounds(self):
        parent = TNode(self.parent_symbol, 0)
        child = TNode(self.child_symbol, 1)
        edge = TEdge(parent, child, self.formula)
        # TODO
