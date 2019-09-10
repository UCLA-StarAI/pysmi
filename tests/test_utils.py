import unittest

from pysmt.shortcuts import Symbol, GE, Or, And, LE, Times, Plus
from pysmt.typing import REAL

from pysmi.utils import *


class TestUtils(unittest.TestCase):
    def test_check_clause(self):
        var_x = Symbol("X", REAL)

        f1 = GE(var_x, Real(1))
        f2 = Or(GE(var_x, Real(1)), GE(Real(0), var_x))
        f3 = And(GE(var_x, Real(1)), GE(Real(0), var_x))

        self.assertTrue(check_clause(f1))
        self.assertTrue(check_clause(f2))
        self.assertFalse(check_clause(f3))

    def test_solve_equation(self):
        pass

    def test_merge_intervals(self):
        intervals = [[1, 2], [3, 4], [2, 3]]
        intervals = merge_intervals(intervals)
        self.assertListEqual(intervals, [[1, 4]])

    def test_domains_to_intervals(self):
        # 1 < x < 2, 3 < 2x < 5, 3 < x < 4.
        var_x = Symbol("x", REAL)

        atom1 = And(LE(Real(1), var_x), LE(var_x, Real(2)))
        atom2 = And(LE(Real(3), Times(Real(2), var_x)), LE(Times(Real(2), var_x), Real(5)))
        atom3 = And(LE(Real(3), var_x), LE(var_x, Real(4)))
        formula = Or([atom1, atom2, atom3])

        interval = domains_to_intervals(formula)
        self.assertListEqual(interval,
                             [[1, float(1.5)], [1.5, 2], [2, 2.5], [3, 4]])

        interval = domains_to_intervals(LE(Real(1), var_x))
        self.assertListEqual(interval, [[1, float('inf')]])

    def test_intervals_to_points(self):
        intervals = [[1, 2], [3, 4]]
        points = intervals_to_points(intervals)
        self.assertListEqual(points, [Real(1), Real(2), Real(3), Real(4)])

    def test_get_critical_points(self):
        var_x = Symbol("x", REAL)

        bounds = [Plus(var_x, Real(1)), Real(2), Real(3)]
        domains = And(LE(var_x, Real(1.5)))
        points = get_critical_points(bounds, domains)
        self.assertListEqual(points, [Real(1)])

    def test_find_edge_critical_points(self):
        sender = Symbol("y", REAL)
        recipient = Symbol("x", REAL)
        sender_domains = And(LE(Real(-1.5), sender), LE(sender, Real(1.5)))
        recipient_domains = And(LE(Real(-1), recipient), LE(recipient, Real(1)))
        edge_formula = Or(LE(recipient + Real(1), sender),
                          LE(sender, recipient - Real(1)))
        critical_points = find_edge_critical_points(
            sender_domains, recipient_domains, edge_formula
        )
        cps = set()
        for cp in critical_points:
            cps.add(float(cp))
        self.assertSetEqual(cps, {-1, -0.5, 0.5, 1})

    def test_get_intervals_intersection(self):
        intervals1 = [[0, 2, 2], [3, 4, 1]]
        intervals2 = [[1.5, 3.5, 5]]
        intervals3 = [[1.5, 2, 7], [3, 5, 9]]
        interval_lists = [intervals1, intervals2, intervals3]
        gold_intersection = [[1.5, 2, 14], [3, 3.5, 15]]
        self.assertListEqual(get_intervals_intersection(interval_lists),
                             gold_intersection)

    def test_initiate_bound(self):
        pass

    def test_get_x_coefficient(self):
        var_x = Symbol("x", REAL)
        term1 = Plus(var_x, Real(1))  # x + 1
        term2 = Plus(Real(1), var_x, )  # 1 + x

        atom1 = Times(Real(2), var_x)  # 2 * x
        atom2 = Times(var_x, Real(2))   # x * 2
        atom3 = Times(Real(3), var_x)  # 2 * x

        term3 = Plus(atom1, Real(1))  # 2 * x + 1
        term4 = Plus(atom2, Real(1))  # x * 2 + 1
        term5 = Plus(Real(1), atom1)  # 1 + 2 * x
        term6 = Plus(Real(1), atom2)  # 1 + x * 2

        term7 = Plus([atom1, atom2, atom3, Real(3)])  # 2 * x + 1 +  x * 2 + 1

        self.assertEqual(get_x_coefficient(term1), Real(1))
        self.assertEqual(get_x_coefficient(term2), Real(1))
        self.assertEqual(get_x_coefficient(term3), Real(2))
        self.assertEqual(get_x_coefficient(term4), Real(2))
        self.assertEqual(get_x_coefficient(term5), Real(2))
        self.assertEqual(get_x_coefficient(term6), Real(2))
        self.assertEqual(get_x_coefficient(term7), Real(7))

        # together test get_coefficients
        self.assertEqual(get_coefficients(term1)[var_x], Real(1))
        self.assertEqual(get_coefficients(term2)[var_x], Real(1))
        self.assertEqual(get_coefficients(term3)[var_x], Real(2))
        self.assertEqual(get_coefficients(term4)[var_x], Real(2))
        self.assertEqual(get_coefficients(term5)[var_x], Real(2))
        self.assertEqual(get_coefficients(term6)[var_x], Real(2))
        self.assertEqual(get_coefficients(term7)[var_x], Real(7))

    def test_get_coefficients(self):
        x = Symbol("x", REAL)
        y = Symbol("y", REAL)
        z = Symbol("z", REAL)
        atom1 = Times(Real(2), x)
        atom2 = Times(Real(3), x)
        atom3 = Times(Real(2), y)
        atom4 = Times(Real(7), z)
        term = Plus([atom1, atom2, atom3, atom4, Real(4), Real(0.1)])

        coefficients = get_coefficients(term)
        self.assertEqual(coefficients[x], Real(5))
        self.assertEqual(coefficients[y], Real(2))
        self.assertEqual(coefficients[z], Real(7))

        term2 = Plus(Real(4), Real(0.1))
        coefficients2 = get_coefficients(term2)
        self.assertEqual(coefficients2, {})

    def test_get_constants(self):
        x = Symbol("x", REAL)
        y = Symbol("y", REAL)
        z = Symbol("z", REAL)
        atom1 = Times(Real(2), x)
        atom2 = Times(Real(3), x)
        atom3 = Times(Real(2), y)
        atom4 = Times(Real(7), z)
        term = Plus([atom1, atom2, atom3, atom4, Real(4), Real(0.1), Real(7.0001)])
        const = get_constants(term)
        const = float(const.constant_value())
        self.assertAlmostEqual(const, 11.1001)

    def test_atom_to_bound(self):
        # univariate atom case
        x = Symbol("x", REAL)
        atom1 = Times(Real(2), x)  # 2x
        atom2 = Times(Real(3), x)  # 3x
        lhs = Plus(atom1, Real(3))
        rhs = Plus(atom2, Real(5))
        atom = LE(lhs, rhs)
        bound, bound_type = atom_to_bound(atom, x)
        self.assertEqual(bound, Real(-2))
        self.assertEqual(bound_type, LOWER)

        x = Symbol("x", REAL)
        atom1 = Times(Real(2), x)
        lhs = Real(3)
        rhs = Plus(atom1, Real(1))
        atom = LE(lhs, rhs)
        bound, bound_type = atom_to_bound(atom, x)
        self.assertEqual(bound, Real(1))
        self.assertEqual(bound_type, LOWER)

        y = Symbol("y", REAL)
        lhs = Plus([y, 2*y, Real(5)])
        rhs = Plus([y, Real(2), Real(9)])
        atom = LE(lhs, rhs)
        bound, bound_type = atom_to_bound(atom, y)
        self.assertEqual(bound, Real(3))
        self.assertEqual(bound_type, UPPER)

        # bivariate case
        lhs = Plus([2*x, 3*y, Real(5)])
        rhs = Plus([3*x, 3*y, Real(4)])
        atom = LE(lhs, rhs)
        bound, bound_type = atom_to_bound(atom, x)
        self.assertEqual(bound, Real(1))
        self.assertEqual(bound_type, LOWER)

        lhs = Plus([2 * x, 3 * y, Real(5)])
        rhs = Plus([3 * x, 4 * y, Real(1)])
        atom = LE(lhs, rhs)
        bound, bound_type = atom_to_bound(atom, x)
        self.assertEqual(bound_type, LOWER)
        y_coef = get_x_coefficient(bound)
        const = get_constants(bound)
        self.assertEqual(y_coef, Real(-1))
        self.assertEqual(const, Real(4))

    def test_get_distinct_vals(self):
        interval = [0, 1]
        self.assertListEqual(get_distinct_vals(interval, 0), [0.5])
        self.assertListEqual(get_distinct_vals(interval, 1), [0.25, 0.5])
        self.assertListEqual(get_distinct_vals(interval, 2), [0.25, 0.5, 0.75])
        self.assertListEqual(get_distinct_vals(interval, 3),
                             [0.125, 0.25, 0.5, 0.75])
        self.assertListEqual(get_distinct_vals(interval, 4),
                             [0.125, 0.25, 0.375, 0.5, 0.75])
        self.assertListEqual(get_distinct_vals(interval, 5),
                             [0.125, 0.25, 0.375, 0.5, 0.625, 0.75])
        self.assertListEqual(get_distinct_vals(interval, 6),
                             [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875])

    def test_interpolate_and_integrate(self):
        xs = [1, 2, 0.5]
        ys = [2, 5, 1.25]
        interval = [0, 3]
        self.assertAlmostEqual(interpolate_and_integrate(xs, ys, interval),
                               12, delta=1e-8)

    def test_path_relabel(self):
        n = 7
        gold_mapping = {3: 0, 1: 1, 5: 2, 0: 3, 2: 4, 4: 5, 6: 6}
        gold_level = {3: 0, 1: 1, 5: 1, 0: 2, 2: 2, 4: 2, 6: 2}
        mapping, level = path_relabel(n)
        self.assertDictEqual(mapping, gold_mapping)
        self.assertDictEqual(level, gold_level)

    def test_half_relabel(self):
        n = 7
        gold_mapping = {0: 6, 1: 4, 2: 2, 3: 0, 4: 1, 5: 3, 6: 5}
        mapping = half_chain_relabel(n)
        self.assertDictEqual(mapping, gold_mapping)

    def test_find_symbolic_bounds(self):
        x = Symbol("x", REAL) # sender
        y = Symbol("y", REAL) # recipient
        sender_constraints = [[0, 1, 1], [2, 3, 1]]  # put 1 instead of p(x)
        test_point = 0.5
        formula = Or(LE(y, x), LE(y, x - 2))
        symbolic_bounds = find_symbolic_bounds(
            sender_constraints, x, test_point, formula
        )
        symbolic_bounds.sort(key=lambda p: p[2], reverse=False)
        sym_bdd1 = symbolic_bounds[0] # should be [y, 1]
        const1 = get_constants(sym_bdd1[0])
        y_coef1 = get_x_coefficient(sym_bdd1[0])
        self.assertEqual(float(const1.constant_value()), 0)
        self.assertEqual(float(y_coef1.constant_value()), 1)
        self.assertEqual(float(sym_bdd1[1].constant_value()), 1)
        sym_bdd2 = symbolic_bounds[1] # should be [2,3]
        self.assertEqual(float(sym_bdd2[0].constant_value()), 2)
        self.assertEqual(float(sym_bdd2[1].constant_value()), 3)





