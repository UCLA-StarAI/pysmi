#Search-based Model Integration (SMI)

Source code for UAI19 paper "Efficient Search-Based Weighted Model Integration". To reproduce the experimental results in the paper, you can run:

```shell script
python3 exapmle.py --graph GRAPH_TYPE --n_nodes NUMBER_OF_NODES [--output_dir DIR]
```    
    
with GRAPH_TYPE being STAR, FAT or PATH which correspond to running model integration on the SMT(lra) theories with star primal graphs, full-3 ary tree primal graphs and path primal graph respectively;
NUMBER_OF_NODES being the maximum number of nodes in the graph to which SMI will run from case with one node;
DIR being optional output directory, which is '.' by default. 
Results and runtime will be output to a .txt file '{GRAPH_TYPE}_{NUMBER_OF_NODES}.txt'.

Dependencies: numpy, scipy, networkx, pysmt.

If you need any help with the code, please contact zhezeng@cs.ucla.edu.


