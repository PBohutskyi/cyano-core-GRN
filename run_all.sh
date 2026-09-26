#!/bin/bash
# run_all.sh: regenerate every table, figure and Dataset S5 to S8 from the inputs.
#
#   cd cyano-core-GRN
#   bash run_all.sh            steps 01 to 14 and Datasets S5 to S8, about two minutes,
#                              from the shipped seeded GENIE3 runs
#   bash run_all.sh --genie3   step 00 first: reruns the 40 seeded GENIE3 runs
#                              (about 30 hours on one core) and stops if any
#                              differs from the shipped run
#
# The four networks of record are fixed inputs and are never regenerated.
# Every step writes logs/<step>.log. The run stops at the first error.

set -e
set -o pipefail
cd "$(dirname "$0")"
mkdir -p outputs/tables outputs/figures logs datasets
PY="${PYTHON:-python3}"

say() { echo "[$(date '+%H:%M:%S')] $*"; }
step() { local name="$1"; shift; say "$name"; "$@" > "logs/$name.log" 2>&1 || { echo "FAILED: $name, see logs/$name.log"; tail -20 "logs/$name.log"; exit 1; }; }

say "software versions"
{ Rscript --version 2>&1
  Rscript -e 'for (p in c("igraph","readxl","openxlsx")) cat(p, as.character(packageVersion(p)), "\n")'
  $PY --version 2>&1
  $PY -c "import pandas, networkx, scipy, matplotlib, openpyxl; print('pandas', pandas.__version__, '| networkx', networkx.__version__, '| scipy', scipy.__version__, '| matplotlib', matplotlib.__version__, '| openpyxl', openpyxl.__version__)"
} > logs/00_versions.log 2>&1
cat logs/00_versions.log

if [ "${1:-}" = "--genie3" ]; then
  [ -d expression ] || { echo "STOP: expression/ is missing"; exit 1; }
  for n in core PCC7942 PCC6803 PCC7002; do
    step 00_genie3_seeded_$n  Rscript scripts/inference/00_genie3_seeded.R $n 1:10
  done
fi

step 01_regulator_counts            Rscript scripts/analysis/01_regulator_counts.R
step 02_figure1A                    Rscript scripts/analysis/02_figure1A.R
step 03_table1_enrichment           Rscript scripts/analysis/03_table1_enrichment.R
step 04_figure2                     Rscript scripts/analysis/04_figure2.R
step 05_centralities                Rscript scripts/analysis/05_centralities.R
step 06_table3_network_properties   Rscript scripts/analysis/06_table3_network_properties.R
step 07_louvain_check               $PY     scripts/analysis/07_louvain_check.py
step 08_figure3B_table              $PY     scripts/analysis/08_figure3B_table.py
step 09_figure3C                    $PY     scripts/analysis/09_figure3C.py
step 10_IC_enrichment_robustness    Rscript scripts/analysis/10_IC_enrichment_robustness.R
step 11_figure4_core_stability      Rscript scripts/analysis/11_figure4_core_stability.R
step 12_figureS2_species_stability  Rscript scripts/analysis/12_figureS2_species_stability.R
step 13_crossvalidation             Rscript scripts/analysis/13_crossvalidation.R
step 14_figures5_6                  $PY     scripts/analysis/14_figures5_6.py
step build_dataset_S5               $PY     scripts/datasets/build_dataset_S5.py
step build_dataset_S6               Rscript scripts/datasets/build_dataset_S6.R
step build_dataset_S7               Rscript scripts/datasets/build_dataset_S7.R
step build_dataset_S8               Rscript scripts/datasets/build_dataset_S8.R
say "done"
