from collections import defaultdict
from pysmt.fnode import FNode
from pysmt.shortcuts import And, Real, simplify, LE, Or
from pysmi.utils import (check_clause, solve_equation,
                         domains_to_intervals, intervals_to_points,
                         get_critical_points, get_real_variables,
                         categorize_bounds, get_bounding_box)


class TNode:
    def __init__(self,
                 symbol: FNode,
                 label: int,
                 domains: FNode = None,
                 cache: bool = True,
                 n_bbox: int = 2,
                 xrange: float = 1):
        """
        :param symbol:
        :param label:
        :param domains:
        :param cache:
        :param n_bbox: the number of disjuctive intervals in domain
        :param xrange:  the range of the whole bounding box,
        from -xrange to xrange
        """
        self.symbol = symbol
        self.label = label
        self.edges = defaultdict(list)  # key label to value TEdge
        self.children = []  # a list of TNode
        self.domains = domains or And(symbol >= Real(-1), symbol <= Real(1))
        if domains:
            self.domains = domains
        else:
            intervals = get_bounding_box(n_bbox, xrange)
            clauses = []
            pre_start, pre_end = intervals[0]
            clauses.append(LE(Real(pre_start), symbol))
            for start, end in intervals[1:]:
                clauses.append(Or(LE(symbol, Real(pre_end)),
                                  LE(Real(start), symbol)))
                pre_end = end
            clauses.append(LE(symbol, Real(pre_end)))
            self.domains = And([c for c in clauses])

        self.initiations = []
        if cache:
            self.cache_index = []
            self.cache = dict()
            self.cache_ide = []  # cache intervals and degree

    def add_edge(self, edge):
        labels = {child.label: child for child in self.children}
        child = edge.child
        if edge.child.label not in labels:
            self.children.append(child)
        self.edges[child.label].append(edge)

    def get_value(self, assignment):
        cache_index = self.cache_index
        cache_index.sort(reverse=False)
        value_key = [assignment[i] for i in cache_index]
        value_key = tuple(value_key)
        if value_key in self.cache:
            return self.cache[value_key]
        else:
            return False

    def update_domains(self, domains):
        domains = simplify(domains)
        self.domains = domains

    def store_value(self, assignment, value):
        cache_index = self.cache_index
        cache_index.sort(reverse=False)
        value_key = [assignment[i] for i in cache_index]
        value_key = tuple(value_key)
        self.cache[value_key] = value

    def add_initiations(self, initiations):
        self.initiations = self.initiations + initiations
        self.initiations = list(set(self.initiations))


class TEdge:
    def __init__(self, parent, child, formula):
        """
        :param formula: clause must be of form that child variable on one side
        and other variables and constant on the other side, e.g. root+c <= child
        :param child: TNode object
        """
        self.parent = parent
        self.child = child
        self.formula = simplify(formula)
        self.check_formula()
        self.upper_bounds, self.lower_bounds, self.critical_points = self.collect_bounds()

    def check_formula(self):
        """
        check if formula is valid for edge representation
        :return:
        """
        assert check_clause(self.formula)
        real_variables = list(get_real_variables(self.formula))
        assert len(real_variables) == 2
        if len(self.formula.get_atoms()) == 1:  # atom
            assert self.formula.is_le() or self.formula.is_lt()
        else:
            for atom in self.formula.args():
                assert atom.is_le() or atom.is_lt()

    def collect_bounds(self):
        """
        collect integration bounds for child variables and critical points
        for root variable
        :return: lists of integration upper bounds, integration lower bounds,
        and critical points, all with type FNode
        """
        upper_bounds, lower_bounds, parent_atoms = categorize_bounds(
            self.formula, self.child.symbol
        )
        critical_points = []

        # collect critical points from parent atoms, i.e. atoms with parent only
        for atom in parent_atoms:
            lower, upper = atom.arg(0), atom.arg(1)
            solution = solve_equation(lower, upper, domains=self.parent.domains)
            if solution:
                critical_points.append(solution[self.parent.symbol])

        # collect critical points from parent domain
        intervals = domains_to_intervals(self.parent.domains)
        critical_points.extend(intervals_to_points(intervals))

        # collect critical points from solving bounds for child
        intervals = domains_to_intervals(self.child.domains)
        points = intervals_to_points(intervals)
        bounds = list(set(upper_bounds + lower_bounds + points))
        critical_points.extend(
            get_critical_points(bounds, domains=self.parent.domains))
        critical_points = list(set(critical_points))

        return upper_bounds, lower_bounds, critical_points
