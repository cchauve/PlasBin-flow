#!/bin/bash
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --account=def-chauvec
#SBATCH --output=/home/chauvec/projects/ctb-chauvec/PLASMIDS/USRA_PLASMIDS_2023/bin/PlasBin-flow/example/tuning/tuning.out
#SBATCH --error=/home/chauvec/projects/ctb-chauvec/PLASMIDS/USRA_PLASMIDS_2023/bin/PlasBin-flow/example/tuning/tuning.err
#SBATCH --job-name=plasbin_tuning

module load StdEnv/2020  gcc/9.3.0 blast+/2.12.0
module load python/3

rm -rf output
mkdir -p output

python ../../code/plasbin_utils.py pls_genes_db \
       --input_file input.csv \
       --out_dir output \
       --tmp_dir tuning_tmp

python ../../code/plasbin_utils.py tuning \
       --input_file    input.csv \
       --out_dir       output \
       --tmp_dir       tuning_tmp \
       --log_file      tuning.log \
       --db_file       output/pls.genes.fasta \
       --out_file      tuning_output.csv \
       --pid_threshold 0.9 \
       --cov_threshold 0.7 \
       --n_gcints      5
