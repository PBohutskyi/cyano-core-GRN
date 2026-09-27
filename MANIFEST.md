# Manifest

Scripts in run order (`run_all.sh`). Paths are relative to the repository root. `tables/` and `figures/` mean `outputs/tables/` and `outputs/figures/`.

## Analysis

| Script | Reads | Writes | Supports |
|---|---|---|---|
| `01_regulator_counts.R` | Datasets S3, S4; `curation/S5_curated_regulators.xlsx` | `tables/ortholog_groups.csv`, `figure2_pipeline_counts.csv`, `figure2_tf_venn.csv` | Figure 2 counts |
| `02_figure1A.R` | Dataset S3; `reference/genome_gene_lists.csv` | `tables/figure1A_regions.csv`; `figures/fig1_panelA` | Figure 1A |
| `03_table1_enrichment.R` | Dataset S3; `reference/ref_protein_coding_7942.csv` | `tables/table1_enrichment.csv`, `table1_significant.csv` | Table 1 |
| `04_figure2.R` | output of 01 | `figures/fig2_*` | Figure 2 |
| `05_centralities.R` | `networks/*/network_of_record.graphml`, `networks/*/seeded/`; curation | `tables/centralities_{core,s7942,s6803,s7002}.csv` | IC rankings; Dataset S7 |
| `06_table3_network_properties.R` | networks of record; regulator lists; Dataset S3; gene lists | `tables/table3_network_properties.csv` | Table 3; core GRN properties in Results |
| `07_louvain_check.py` | core network of record | `tables/louvain_modules.csv`, `louvain_agreement.csv`, `louvain_run_log.txt` | Louvain modules of Figure 3A |
| `08_community_enrichment.R` | core network of record (Louvain communities); `reference/ref_protein_coding_7942.csv`; curation | `tables/community_enrichment.csv`, `community_enrichment_significant.csv` | COG and KEGG enrichment of the communities: Figure 3B, Results 3.4.1 |
| `09_figure3B_table.py` | `tables/centralities_core.csv`, `community_enrichment_significant.csv`; curation | `tables/figure3B_table.xlsx` | Figure 3B |
| `10_figure3C.py` | core network of record; `tables/centralities_core.csv`; curation | `figures/fig3C_assembled`; `tables/fig3C_values.csv`, `fig3C_leaf_order.txt` | Figure 3C |
| `11_IC_enrichment_robustness.R` | `tables/centralities_core.csv`; curation | `tables/S6-4_IC_robustness.csv`, `S6-5_IC_enrichment.csv`, `core_top20_IC.csv` | Top-15 enrichment; Dataset S6-4, S6-5 |
| `12_figure4_core_stability.R` | core network of record and seeded runs; `tables/centralities_core.csv` | `figures/fig4_*`; `tables/fig4_*.csv`, `seed_centralities_core.csv`, `seed_outdegree_core.csv`, `fig4_numbers_for_SI.txt` | Figure 4; Dataset S8 |
| `13_figureS2_species_stability.R` | species networks of record and seeded runs; centrality tables; curated S5-4 | `figures/figS2_*`; `tables/figS2_*.csv`, `seed_centralities_species.csv`, `seed_outdegree_species.csv`, `figS2_numbers_for_SI.txt` | Figure S2; Dataset S8 |
| `14_network_recurrence.R` | networks of record and seeded runs | `tables/edge_recurrence.csv`, `network_recurrence_summary.csv` | Recurrence of each network in its seeded runs; Dataset S8-6 |
| `15_crossvalidation.R` | centrality tables; curated S5-4; curation | `tables/S6-6_xval_enrichment.csv`, `S6-7_xval_sensitivity.csv`, `xval_instances.csv` | Cross-validation; Dataset S6-6 to S6-8 |
| `16_figures5_6.py` | centrality tables; `tables/ortholog_groups.csv`; curation | `tables/figure5_symbols.csv`, `figure6_symbols.csv`, `figure6_regulators.csv`, `figure56_counts.txt` | Figures 5 and 6 (placements) |

## Datasets

| Script | Writes |
|---|---|
| `build_dataset_S5.py` | `datasets/Dataset_S5_Curated_TFs.xlsx` |
| `build_dataset_S6.R` | `datasets/Dataset_S6_Core_GRN_nodes_TFs_and_edges.xlsx` |
| `build_dataset_S7.R` | `datasets/Dataset_S7_Centrality_and_stability.xlsx` |
| `build_dataset_S8.R` | `datasets/Dataset_S8_Seeded_run_stability.xlsx` |

## Step 00, seeded GENIE3 (`bash run_all.sh --genie3`)

| Script | Reads | Writes |
|---|---|---|
| `00_genie3_seeded.R` | `expression/expression_<network>.csv`; `networks/<network>/regulators_GENIE3_input.csv` | `outputs/genie3_rerun/<network>/`, with a comparison against `networks/<network>/seeded/` |

## Inputs

| File | Contents |
|---|---|
| `networks/<network>/network_of_record.graphml` | Network of record: core GRN 889 nodes, 1,100 edges; PCC 7942 1,978 and 3,500; PCC 6803 2,461 and 3,500; PCC 7002 2,181 and 3,500 |
| `networks/<network>/regulators_GENIE3_input.csv` | Regulators supplied to GENIE3: core 42, PCC 7942 71, PCC 6803 74, PCC 7002 78 |
| `networks/<network>/seeded/top<N>_seed01..10.csv` | Top edges of the ten seeded GENIE3 runs (regulatoryGene, targetGene, weight) |
| `networks/<network>/seeded/run_manifest.csv`, `run_timing.csv` | Settings, versions, input size and run time of the seeded runs |
| `expression/expression_<network>.csv` | Log2 TPM expression compendia supplied to GENIE3 (genes in rows, samples in columns) |
| `curation/core_regulator_curation.csv` | Names and functional classification of the 38 core regulators |
| `curation/S5_curated_regulators.xlsx` | Curated regulator lists per strain; conserved regulators and orthologs (S5-4) |
| `curation/S6_source_annotation.xlsx` | Gene, regulator and edge annotation of the core GRN |
| `curation/figure6_categories.csv` | Names and functional categories of the Figure 6 regulators |
| `reference/ref_protein_coding_7942.csv` | PCC 7942 protein-coding genes with COG and KEGG annotation |
| `reference/genome_gene_lists.csv` | Protein-coding genes of the three genomes |
