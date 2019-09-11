import os

from pysmi.graph2formula import *
from pysmi.smi import SMI
import time
from tqdm import tqdm
import argparse


def run_example(graph_type,
                n_nodes,
                output_dir):
    output_path = os.path.join(
        output_dir,
        f"{graph_type}_{n_nodes}.txt"
    )
    f = open(output_path, 'w')
    runtimes = []
    results = []
    for i in tqdm(range(1, n_nodes + 1)):
        graph = create_graph(graph_type, i)
        stree = graph_to_tree(graph)

        t1 = time.time()
        result = SMI.compute_mi(stree)
        t2 = time.time()
        results.append(result)
        runtimes.append((t2 - t1))
        f.write("with number of node {}, model "
                "integration result: {}\n".format(i, result))
        f.write("smi runtime: {}\n".format(t2 - t1))
        f.write("results so far: {}\n".format(results))
        f.write("runtime so far: {}\n".format(runtimes))
        f.flush()
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--n_nodes", type=int)
    parser.add_argument("-o", "--output_dir", default='.')
    parser.add_argument(
        "-g", "--graph",
        help="primal graph, STAR, FAT, or PATH, CHAIN, HALF_CHAIN, STREET"
    )
    args = parser.parse_args()
    run_example(
        args.graph,
        args.n_nodes,
        args.output_dir
    )
