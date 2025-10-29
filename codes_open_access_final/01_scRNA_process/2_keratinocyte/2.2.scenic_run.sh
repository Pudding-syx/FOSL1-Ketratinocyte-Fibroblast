#!/bin/bash
#SBATCH -J pyscenic
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --exclusive
#SBATCH -n 48
#SBATCH --mem 256g
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err



pyscenic grn --num_workers 12 --output adj.sample.tsv --method grnboost2 2.2.scenic_input.loom allTFs_hg38.txt
pyscenic ctx adj.sample.tsv hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname 10.4.scenic_input.loom --mode "dask_multiprocessing" --output reg.csv --num_workers 12 --mask_dropouts
pyscenic aucell 2.2.scenic_input.loom reg.csv --output 2.2.scenic_result.loom --num_workers 12
