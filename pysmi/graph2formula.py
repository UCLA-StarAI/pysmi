import networkx as nx
from networkx import Graph
from pysmt.fnode import FNode
from pysmt.shortcuts import Symbol, Or, Real, And, LE
from pysmt.typing import REAL
from pysmi.smtree import TNode, TEdge
from pysmi.utils import domains_to_intervals, is_literal, half_chain_relabel


def create_graph(graph_type, n_nodes, n_branches=3):
    """
    generate graph from networkx
    :param graph_type: STAR, FAT/SNOW, PATH, CHAIN
    :param n_nodes: number of nodes
    :param n_branches: number of branches in snow case
    :return:
    """
    if graph_type == 'PATH':
        graph = nx.generators.path_graph(n_nodes)
    elif graph_type == 'FAT' or graph_type == 'SNOW':
        graph = nx.full_rary_tree(n_branches, n_nodes)
    elif graph_type == 'STAR':
        graph = nx.star_graph(n_nodes)
    elif graph_type == 'CHAIN':
        graph = nx.generators.path_graph(n_nodes)
        mapping = half_chain_relabel(n_nodes)
        graph = nx.relabel_nodes(graph, mapping=mapping)
    else:
        graph = None
    return graph


def graph_to_tree(graph: Graph,
                  constant: float = 1.0,
                  return_map=False,
                  n_bbox=1):
    edges, nodes = list(graph.edges), list(graph.nodes)
    n = len(nodes)
    node_map = {}
    for label in range(n):
        symbol = Symbol("x_{}".format(label), REAL)
        if n_bbox:
            node_map[label] = TNode(label=label, symbol=symbol, n_bbox=n_bbox)
        else:
            node_map[label] = TNode(label=label, symbol=symbol)
    for u, v in edges:
        u, v = min(u, v), max(u, v)  # u -> v
        parent = node_map[u]
        child = node_map[v]
        formula = Or(parent.symbol + Real(constant) <= child.symbol,
                     parent.symbol + Real(-constant) >= child.symbol)
        edge = TEdge(parent=parent, child=child, formula=formula)
        parent.add_edge(edge)

    if return_map:
        return node_map[0], node_map

    return node_map[0]


def tree_to_formula(tree: TNode) -> FNode:
    """
    output SMT formula represented by the input SMT tree
    :param tree: SMT tree
    :return: SMT formula in CNF form
    """
    from collections import deque
    formulas = []  # list of clauses
    queue = deque([tree])
    while queue:
        node = queue.popleft()
        edges = []
        for child in node.children:
            queue.append(child)
            edges.extend(node.edges[child.label])

        for e in edges:
            f = e.formula
            if f.is_or() or is_literal(f):
                formulas.append(f)
            else:
                for c in f.args():
                    assert is_literal(c) or c.is_or()
                    formulas.append(c)

        if node.domains:
            intervals = domains_to_intervals(node.domains)
            clauses = []
            # if intervals:
            #     # print('NNNN', node.domains)
            pre_start, pre_end = intervals[0]
            clauses.append(LE(Real(pre_start), node.symbol))
            for start, end in intervals[1:]:
                clauses.append(Or(LE(node.symbol, Real(pre_end)),
                                  LE(Real(start), node.symbol)))
                pre_end = end
            clauses.append(LE(node.symbol, Real(pre_end)))
            formulas.extend(clauses)
    return And([c for c in formulas]) if formulas else None


def check_tree(root: TNode):
    """
    print node and edge information in the input SMT tree
    :param root: SMT tree
    :return: none
    """
    from collections import deque
    que = deque([root])
    checked = []
    print("-----------------------------------------------")
    while que:
        node = que.popleft()
        # queue = queue + node.child
        checked.append(node.label)
        for child in node.children:
            if child.label not in checked:  # in case loop
                que.append(child)
        print("tree root has symbol {}".format(node.symbol))
        print("tree root has label {}".format(node.label))
        print("tree root has domains {}".format(node.domains))
        if not node.children:
            print("tree with root {} has no child".format(node.symbol))
        else:
            for child in node.children:
                print("-----------------------------------------------")
                edges = node.edges[child.label]
                print("child:{}".format(child.symbol))
                i = 0
                if not edges:
                    print("no edge between root {} and child {}".format(
                        node.symbol, child.symbol))
                else:
                    for edge in edges:
                        i = i + 1
                        print("edge {}:".format(i))
                        print("edge represents formula {}".format(edge.formula))
                        print("edge with child {}".format(edge.child.symbol))
                        if edge.critical_points:
                            print("critcal points {}".format(
                                edge.critical_points
                            ))
                        if edge.upper_bounds:
                            print("upper bounds {}".format(edge.upper_bounds))
                        if edge.lower_bounds:
                            print("lower bounds {}".format(edge.lower_bounds))
            print("-----------------------------------------------")


def test_tree():
    r = 3
    n = 8
    graph = nx.generators.full_rary_tree(r, n)
    root = graph_to_tree(graph)
    check_tree(root)


if __name__ == '__main__':
    test_tree()
