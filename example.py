import networkx as nx

from pysmi.graph2formula import *
from pysmi.utils import half_chain_relabel
from pysmi.smi import SMI

from pysmt.shortcuts import *
from pysmt.typing import REAL

import time
from tqdm import tqdm
import argparse

if __name__ == '__main__':
    graph_type = "STAR"
    n_nodes = 6
    jump = n_nodes - 1
    cache = True

    for i in tqdm(range(jump, n_nodes)):
        if graph_type == 'PATH':
            graph = nx.generators.path_graph(i)
        elif graph_type == 'FAT':
            n_branches = 3  # default
            graph = nx.full_rary_tree(n_branches, i)
        elif graph_type == 'STAR':
            graph = nx.star_graph(i)
        elif graph_type == 'CHAIN':
            graph = nx.generators.path_graph(i)
            mapping = half_chain_relabel(i)
            graph = nx.relabel_nodes(graph, mapping=mapping)
        else:
            print("not valid graph type!")
            break

        stree = graph_to_tree(graph)
        formula = tree_to_formula(stree)
        check_tree(stree)

        e1 = time.time()
        result = SMI.compute_mi(stree)
        e2 = time.time()

        print("mi result: {}".format(result))
        print("time: {}".format(e2 - e1))
