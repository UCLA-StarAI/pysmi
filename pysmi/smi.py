from collections import defaultdict

from pysmi.graph2formula import tree_to_formula
from pysmi.smtree import TNode, TEdge
from pysmi.utils import *


class SmiTask:
    def __init__(self, node: TNode, cache):
        self.node = node
        self.formula = tree_to_formula(self.node)
        self.intervals = domains_to_intervals(self.node.domains)
        self.progress = 0
        self.cache = cache
        self.res_map = defaultdict(dict)
        self.reversed_info = None
        self.sub_task_info = self._generate_sub_task_info()

    def _generate_sub_task_info(self):
        # to get intervals and degrees
        interval_and_degree = [[s, t, 0] for s, t in self.intervals]
        interval_and_degree = get_intervals_intersection(
            [interval_and_degree, self.node.cache_ide]
        )
        label = self.node.label
        sub_tasks = []
        for s, t, d in interval_and_degree:
            # for x in get_distinct_vals([s, t], d):
            for x in get_distinct_vals_cache([s, t], d, self.cache, label):
                # If the a specific value of x for a particular node
                # is cached, then no need to calculate its children
                if (label, x) in self.cache:
                    self.res_map[(s, t)][x] = self.cache[(label, x)]
                    continue
                for child in self.node.children:
                    sub_tasks.append([s, t, d, x, child])
        return sub_tasks

    def has_sub_tasks(self):
        return self.progress < len(self.sub_task_info)

    def add_result(self, res):
        s, t, _, x, child = self.sub_task_info[self.progress]
        self.res_map[(s, t)][x] = self.res_map[(s, t)].get(x, 1) * res

    def generate_next_task(self):
        s, t, d, x, child = self.sub_task_info[self.progress]
        self.reversed_info = child.domains
        for edge in self.node.edges[child.label]:
            edge_formula = simplify(substitute(
                edge.formula,
                {self.node.symbol: Real(x)})
            )
            domains = child.domains.And(edge_formula)
            child.update_domains(domains)
        return SmiTask(child, self.cache)

    def finish_current_task(self):
        _, _, _, _, child = self.sub_task_info[self.progress]
        child.domains = self.reversed_info
        self.reversed_info = None
        self.progress += 1

    def compute_1d_formula(self):
        assert not get_boolean_variables(self.formula) and \
            len(get_real_variables(self.formula)) == 1
        result = [abs(end - start) for start, end in self.intervals]
        return sum(result)

    def summarize(self):
        if not self.node.children:
            return self.compute_1d_formula()
        sum_ = 0
        label = self.node.label
        for (s, t), vals in self.res_map.items():
            xs, ys = [], []
            for x, y in vals.items():
                self.cache[(label, x)] = y
                xs.append(x)
                ys.append(y)
            sum_ += interpolate_and_integrate(xs, ys, [s, t])
        return sum_


class SMI:

    @staticmethod
    def compute_mi(tree):
        SMI.pe_node(tree)
        return SMI.compute_mi_aos(tree)

    @staticmethod
    def compute_mi_aos(tree):
        cache = {}
        root_task = SmiTask(tree, cache)
        task_stack = [root_task]
        while True:
            task = task_stack[-1]
            if not task.has_sub_tasks():
                result = task.summarize()
                if task.node.label == 0:
                    return result
                task_stack.pop()
                parent = task_stack[-1]
                parent.add_result(result)
                parent.finish_current_task()
            else:
                task_stack.append(task.generate_next_task())

    @staticmethod
    def pe_edge(edge: TEdge,
                child_constraints=None):
        """
        obtain integration bounds and degree for root variable from edge formula
        and constraints on child variable
        :param edge: TEdge object
        :param child_constraints: constraints for child variable,
        lists of elements of form [start pint, end point, integration degree], type float
        :return: integration intervals and degree,
        lists of elements of form [start pint, end point, integration degree], type float
        """
        upper_bounds, lower_bounds = edge.upper_bounds, edge.lower_bounds
        critical_points = edge.critical_points

        formula = edge.formula.And(edge.parent.domains).And(edge.child.domains)
        if not child_constraints:
            child_domains_bounds = domains_to_intervals(edge.child.domains)
            domain_bounds = []
            for bound in child_domains_bounds:
                domain_bounds.append(
                    [bound[0], bound[1], 0])  # set default degree to 0
        else:
            domain_bounds = child_constraints
            child_symbol = edge.child.symbol
            bounds = upper_bounds + lower_bounds
            child_domains = []
            for constraint in child_constraints:
                start, end = constraint[0], constraint[1]
                child_domains.append(
                    And(child_symbol > Real(start), child_symbol < Real(end)))
                bounds = bounds + [Real(start), Real(end)]
            child_domains = Or([atom for atom in child_domains])
            formula = formula.And(child_domains)
            critical_points = critical_points + get_critical_points(bounds,
                                                                    domains=edge.parent.domains)

        # append critical points from root domain
        root_domain_bounds = domains_to_intervals(edge.parent.domains)
        root_domain_points = intervals_to_points(root_domain_bounds)
        critical_points = critical_points + root_domain_points

        critical_points = list(set(critical_points))
        formula = simplify(formula)
        critical_points = [float(point.constant_value()) for point in
                           critical_points]
        critical_points.sort(reverse=False)

        interval_and_degree = []
        # not include inf, for bounded variables only
        intervals = zip(critical_points, critical_points[1:])
        root_symbol = edge.parent.symbol
        for start, end in intervals:
            # use mid point as default test point
            test_point = (start + end) * 0.5
            test_interval = Equals(root_symbol, Real(test_point))

            with Solver(name="msat") as solver:
                solver.add_assertion(formula.And(test_interval))
                if not solver.solve():
                    continue

            # test_point = (start + end) * 0.5
            # use mid point as default test point
            num_upper_bounds, num_lower_bounds = dict(), dict()
            for upper in upper_bounds:
                num_upper_bounds[initiate_bound(upper, test_point)] = upper
            for lower in lower_bounds:
                num_lower_bounds[initiate_bound(lower, test_point)] = lower

            max_upper = float('-inf') if not upper_bounds else max(
                list(num_upper_bounds.keys()))
            min_lower = float('inf') if not lower_bounds else min(
                list(num_lower_bounds.keys()))

            integration_bounds = []
            for domain_bound in domain_bounds:
                min_domain, max_domain, degree = domain_bound

                if max_upper >= max_domain or min_lower <= min_domain or max_upper >= min_lower:
                    integration_bounds.append(
                        [Real(min_domain), Real(max_domain), degree])
                else:
                    if min_lower < max_domain:
                        integration_bounds.append(
                            [num_lower_bounds[min_lower], Real(max_domain),
                             degree])
                    if min_domain < max_upper:
                        integration_bounds.append(
                            [Real(min_domain), num_upper_bounds[max_upper],
                             degree])

            degrees = [get_degree(bound[0], bound[1], bound[2]) for bound in
                       integration_bounds]
            interval_and_degree.append([start, end, max(degrees)])
        return interval_and_degree

    @staticmethod
    def pe_node(tree: TNode):
        """
        obtain integration intervals and polynomial degree of
        integration w.r.t. root variable
        :param tree: SMT tree
        :return: integration intervals and degree
        """

        if not tree.children:
            intervals = domains_to_intervals(tree.domains)
            interval_and_degree = []
            for interval in intervals:
                interval.append(0)
                interval_and_degree.append(interval)  # append default degree 0
            return interval_and_degree

        interval_and_degree = []
        for child in tree.children:
            child_constraints = SMI.pe_node(child)
            for edge in tree.edges[child.label]:
                interval = SMI.pe_edge(edge, child_constraints)
                if interval:
                    interval_and_degree.append(interval)

        interval_and_degree = get_intervals_intersection(interval_and_degree)
        tree.cache_ide = interval_and_degree
        return interval_and_degree
