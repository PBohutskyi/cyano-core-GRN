# Conserved core gene regulatory network of three cyanobacteria

Code and data for:

> Bohutskyi P, DiMura R, Johnson Z, Li R, Anderson D, Cheung M. Stress-related transcriptional regulators enriched in conserved core GRN for three cyanobacteria: network topology maps the highest-influence nodes as candidate engineering targets. *Microbiology Spectrum* (under review, manuscript Spectrum01958-26).

Archived release: v1.1, https://doi.org/10.5281/zenodo.23002537

The repository reproduces every table, figure and Supplementary Dataset S5 to S8 of the paper from the networks of record and the curated inputs, in about two minutes. It also holds the expression data and the script to rerun the seeded GENIE3 inference.

## Quick start

Requirements (versions used for the paper; see `logs/00_versions.log` after a run):

- R 4.6.0 with igraph 2.3.3, readxl 1.5.0, openxlsx 4.2.8.1 (and GENIE3 1.34.0 from Bioconductor, only to rerun GENIE3)
- Python 3.13 with networkx 3.4.2, pandas 3.0.5, scipy 1.18.1, matplotlib 3.11.1, openpyxl 3.1.5

```
bash run_all.sh
```

The pipeline was also run on Linux with R 4.3.3 and the same package versions: every output was identical except one value, the Figure 4A Wilcoxon P (0.395 instead of 0.397), because R 4.4 changed how `wilcox.test` treats ties.

All outputs go to `outputs/tables`, `outputs/figures` and `datasets/`; each step writes a log to `logs/`. To use a specific Python, run `PYTHON=/path/to/python bash run_all.sh`.

## What is where

| Folder | Contents |
|---|---|
| `networks/<network>/` | For the core GRN and the three species GRNs (PCC 7942, PCC 6803, PCC 7002): the network of record (`network_of_record.graphml`), the regulators supplied to GENIE3, and the ten seeded GENIE3 runs (`seeded/`) |
| `curation/` | Curated inputs: regulator names and functional classification, species regulator lists and ortholog map, gene and edge annotation, Figure 6 names and categories. No script writes these files |
| `reference/` | PCC 7942 protein-coding annotation (Table 1 background) and the protein-coding gene lists of the three genomes (Figure 1A) |
| `datasets/` | Supplementary Datasets S1 to S8. S1 to S4 are fixed data (S3 and S4 are also script inputs); S5 to S8 are built by the scripts |
| `expression/` | Log2 TPM expression compendia supplied to GENIE3 (genes in rows, samples in columns) |
| `scripts/analysis/` | Analysis scripts, numbered in run order |
| `scripts/datasets/` | Builders of Datasets S5 to S8 |
| `scripts/inference/` | Seeded GENIE3 inference |
| `outputs/` | Tables and figures written by the scripts |

`MANIFEST.md` lists every script with its inputs and outputs and the table, figure or Dataset each output supports.

## Networks

**Networks of record.** Each network was inferred once with GENIE3 without a fixed random seed and pruned to its top edges: 1,111 for the core GRN, then its largest connected component (889 nodes, 1,100 edges), and 3,500 for each species GRN. These four fixed networks are the inputs of every downstream calculation and are not regenerated. They are stored as GraphML: directed, RefSeq locus tags as node identifiers, node attribute `TF` (1 for regulators supplied to GENIE3), and for the core GRN the Louvain module (`louvain_cluster`, Figure 3). Edges carry no weights; all measures are computed on the unweighted networks.

The core GRN was inferred from 1,312 genes and 42 regulators, giving 55,062 regulator-to-gene links before pruning (42 x 1,311); `networks/core/seeded/run_manifest.csv` and `run_timing.csv` record these values for the seeded runs.

**Seeded runs.** To test reproducibility (Figure 4, Figure S2, Dataset S8), each network was inferred ten more times with GENIE3 and a fixed seed (1 to 10), with the same input and settings, and pruned in the same way. The top edges of each run are in `networks/<network>/seeded/`.

## Rerunning GENIE3

Step 00 of the pipeline regenerates the 40 seeded runs from the expression data and checks each against the shipped run; the analysis then continues with steps 01 to 16:

```
bash run_all.sh --genie3
```

This takes about 30 hours on one core. One network and seed can also be run on its own:

```
Rscript scripts/inference/00_genie3_seeded.R PCC7002 1
```

The first argument is the network (`core`, `PCC7942`, `PCC6803`, `PCC7002`), the second the seeds (default `1:10`). Results go to `outputs/genie3_rerun/<network>/`, and the script reports whether the top edges are identical to the shipped run. Identical results need GENIE3 1.34.0 and one core. One run takes about 20 minutes (PCC 7002), 30 (PCC 7942), 50 (core) or 80 (PCC 6803) on one core.

## Licence and citation

Code: MIT licence (`LICENSE`). Data (curated inputs, networks, expression compendia and Datasets): CC BY 4.0 (`LICENSE-DATA.md`). Please cite the paper when you use this material; `CITATION.cff` gives the reference.
