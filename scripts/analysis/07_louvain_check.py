#!/usr/bin/env python3
"""
07_louvain_check.py

Checks that the Louvain partition stored in the core GRN (node attribute
louvain_cluster; Figure 3A, Dataset S6-2) is reproduced by rerunning Louvain.

Method: NetworkX 3.4.2, Louvain on the directed core network, resolution 5,
seed 12, modules smaller than 9 nodes discarded. Other NetworkX versions may
return a different partition.

Input  : networks/core/network_of_record.graphml
Outputs: outputs/tables/louvain_modules.csv     one row per node
         outputs/tables/louvain_agreement.csv   rerun against the stored partition
         outputs/tables/louvain_run_log.txt     versions, parameters, agreement
"""

import csv
import os
import sys

import networkx as nx

RESOLUTION = 5
SEED = 12
MIN_MODULE_SIZE = 9

GRAPHML = os.path.join("networks", "core", "network_of_record.graphml")
OUT = os.path.join("outputs", "tables")


def adjusted_rand(a, b):
    """Adjusted Rand index between two labellings. No sklearn dependency."""
    from collections import Counter
    from math import comb

    pairs = Counter(zip(a, b))
    ra = Counter(a)
    rb = Counter(b)
    n = len(a)
    sum_ij = sum(comb(v, 2) for v in pairs.values())
    sum_a = sum(comb(v, 2) for v in ra.values())
    sum_b = sum(comb(v, 2) for v in rb.values())
    expected = sum_a * sum_b / comb(n, 2)
    maximum = 0.5 * (sum_a + sum_b)
    return (sum_ij - expected) / (maximum - expected)


def main():
    if nx.__version__ != "3.4.2":
        print("WARNING: NetworkX %s is installed, not 3.4.2. The partition of "
              "record was produced with 3.4.2 and other versions may differ."
              % nx.__version__, file=sys.stderr)

    if not os.path.exists(GRAPHML):
        sys.exit("Input not found: %s" % GRAPHML)
    os.makedirs(OUT, exist_ok=True)

    G = nx.read_graphml(GRAPHML)
    print("read %s: %d nodes, %d edges, directed=%s"
          % (os.path.basename(GRAPHML), G.number_of_nodes(),
             G.number_of_edges(), G.is_directed()))

    communities = nx.community.louvain_communities(
        G, resolution=RESOLUTION, seed=SEED)
    kept = [c for c in communities if len(c) >= MIN_MODULE_SIZE]
    kept.sort(key=len, reverse=True)
    print("%d communities, %d with at least %d nodes"
          % (len(communities), len(kept), MIN_MODULE_SIZE))

    module = {}
    for i, c in enumerate(kept):
        for node in c:
            module[node] = i

    with open(os.path.join(OUT, "louvain_modules.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["locus_tag", "module", "module_size", "is_regulator"])
        for node, data in G.nodes(data=True):
            m = module.get(node, "")
            size = len(kept[m]) if m != "" else ""
            w.writerow([node, m, size, str(data.get("TF", ""))])

    # agreement with the partition stored in the GraphML (Figure 3, Dataset S6-2)
    stored = {n: d["louvain_cluster"] for n, d in G.nodes(data=True)
                 if "louvain_cluster" in d}
    shared = [n for n in stored if n in module]
    ari = adjusted_rand([str(stored[n]) for n in shared],
                        [str(module[n]) for n in shared]) if shared else float("nan")

    with open(os.path.join(OUT, "louvain_agreement.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["locus_tag", "module_rerun", "cluster_of_record", "is_regulator"])
        for n in sorted(shared):
            w.writerow([n, module[n], stored[n],
                        str(G.nodes[n].get("TF", ""))])

    log = os.path.join(OUT, "louvain_run_log.txt")
    with open(log, "w") as fh:
        fh.write("Louvain clustering of the core GRN\n")
        fh.write("python: %s\n" % sys.version.split()[0])
        fh.write("networkx: %s\n" % nx.__version__)
        fh.write("input: %s\n" % GRAPHML)
        fh.write("nodes: %d, edges: %d, directed: %s\n"
                 % (G.number_of_nodes(), G.number_of_edges(), G.is_directed()))
        fh.write("resolution: %s, seed: %s, minimum module size: %s\n"
                 % (RESOLUTION, SEED, MIN_MODULE_SIZE))
        fh.write("communities: %d, retained: %d, nodes covered: %d\n"
                 % (len(communities), len(kept), len(module)))
        fh.write("adjusted Rand index vs the assignment of record: %.4f\n" % ari)
        fh.write("Module integers are assigned in discovery order and carry no "
                 "meaning. An index of 1.00 means the partition is identical "
                 "and only the labels differ.\n")

    print("adjusted Rand index vs the assignment of record: %.4f" % ari)
    print("wrote %s" % os.path.join(OUT, "louvain_modules.csv"))
    print("wrote %s" % os.path.join(OUT, "louvain_agreement.csv"))
    print("wrote %s" % log)


if __name__ == "__main__":
    main()
