#!/bin/bash
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --account=def-chauvec
#SBATCH --output=/home/chauvec/projects/ctb-chauvec/PLASMIDS/USRA_PLASMIDS_2023/bin/PlasBin-flow/code/plasbin_utils_test_noerror.out
#SBATCH --error=/home/chauvec/projects/ctb-chauvec/PLASMIDS/USRA_PLASMIDS_2023/bin/PlasBin-flow/code/plasbin_utils_test_noerror.err
#SBATCH --job-name=plasbin_utils_test_noerror

module load StdEnv/2020  gcc/9.3.0 blast+/2.12.0
module load python/3

INPUT_PREF=./plasbin_utils_test_noerror
TEST_DIR=/scratch/chauvec/PlasBin-flow/plasbin_utils_test_noerror

rm -rf ${TEST_DIR}
mkdir -p ${TEST_DIR}

echo "----> Testing creating plasmids genes database"
python plasbin_utils.py pls_genes_db \
       --input_file ${INPUT_PREF}.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp \
       --log_file   ${TEST_DIR}/pls_genes_db.log

echo "----> Testing mapping plasmid genes to contigs"
python plasbin_utils.py map_genes_to_ctgs \
       --input_file ${INPUT_PREF}.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp \
       --db_file    ${TEST_DIR}/plasbin_utils_test_results/pls.genes.fasta \
       --out_file   ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c.csv \
       --log_file   ${TEST_DIR}/map_genes_to_ctgs.log

echo "----> Testing computing ground truth"
python plasbin_utils.py ground_truth \
       --input_file    ${INPUT_PREF}.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_gt.csv \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --log_file      ${TEST_DIR}/ground_truth_1.log

echo "----> Testing computing gene density with mappings"
python plasbin_utils.py gene_density \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd.csv \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --log_file      ${TEST_DIR}/gene_density_1.log

echo "----> Testing computing GC content intervals"
python plasbin_utils.py gc_intervals \
       --input_file ${INPUT_PREF}.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp \
       --n_gcints   5 \
       --log_file   ${TEST_DIR}/gc_intervals.log

echo "----> Testing computing seeds parameters with gene density and mapping (not needed)"
echo "------> Computing ground truth with gene density and mapping (not needed)"
python plasbin_utils.py ground_truth \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd_gt.csv \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --log_file      ${TEST_DIR}/ground_truth_seeds.log
python plasbin_utils.py seeds \
       --input_file ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd_gt.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp \
       --log_file   ${TEST_DIR}/seeds_1.log

echo "----> Testing tuning with gene density and ground truth"
python plasbin_utils.py tuning \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd_gt.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results_tuning \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp_tuning \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd_gt_tuning.csv \
       --n_gcints      5 \
       --log_file      ${TEST_DIR}/tuning_1.log

echo "----> Testing tuning with gene density and no ground truth"
python plasbin_utils.py tuning \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results_tuning \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp_tuning \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd_tuning.csv \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --n_gcints      5 \
       --log_file      ${TEST_DIR}/tuning_2.log

echo "----> Testing tuning without gene density and ground truth, with mappings"
echo "------> Computing ground truth with mappings"
python plasbin_utils.py ground_truth \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gt.csv \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --log_file      ${TEST_DIR}/ground_truth_tuning.log
python plasbin_utils.py tuning \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gt.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results_tuning \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp_tuning \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gt_tuning.csv \
       --n_gcints      5 \
       --log_file      ${TEST_DIR}/tuning_3.log

echo "----> Testing tuning with without gene density and mappings, with ground truth"
python plasbin_utils.py tuning \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_gt.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results_tuning \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp_tuning \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_gt_tuning.csv \
       --db_file       ${TEST_DIR}/plasbin_utils_test_results/pls.genes.fasta \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --n_gcints      5 \
       --log_file      ${TEST_DIR}/tuning_4.log

echo "----> Testing tuning without gene density and ground truth, with mappings"
python plasbin_utils.py tuning \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results_tuning \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp_tuning \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_tuning.csv \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --n_gcints      5 \
       --log_file      ${TEST_DIR}/tuning_5.log

echo "----> Testing tuning without gene density, ground truth and mappings"
python plasbin_utils.py tuning \
       --input_file    ${INPUT_PREF}.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results_tuning \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp_tuning \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_tuning.csv \
       --db_file       ${TEST_DIR}/plasbin_utils_test_results/pls.genes.fasta \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --n_gcints      5 \
       --log_file      ${TEST_DIR}/tuning_6.log

echo "----> Testing computing GC content probabilities with default intervals"
python plasbin_utils.py gc_probabilities \
       --input_file   ${INPUT_PREF}.csv \
       --out_dir      ${TEST_DIR}/plasbin_utils_test_results_gcp_default \
       --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp \
       --log_file     ${TEST_DIR}/gc_probabilities_0.log

echo "----> Testing computing GC content probabilities"
python plasbin_utils.py gc_probabilities \
       --input_file   ${INPUT_PREF}.csv \
       --out_dir      ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp \
       --gc_intervals ${TEST_DIR}/plasbin_utils_test_results/gc.txt \
       --out_file     ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_gcp.csv \
       --log_file     ${TEST_DIR}/gc_probabilities_1.log

echo "----> Testing computing GC content probabilities with mappings"
python plasbin_utils.py gc_probabilities \
       --input_file   ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c.csv \
       --out_dir      ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp \
       --gc_intervals ${TEST_DIR}/plasbin_utils_test_results/gc.txt \
       --out_file     ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gcp.csv \
       --log_file     ${TEST_DIR}/gc_probabilities_2.log

echo "----> Testing computing GC content probabilities with mappings and gene density"
python plasbin_utils.py gc_probabilities \
       --input_file   ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd.csv \
       --out_dir      ${TEST_DIR}/plasbin_utils_test_results \
       --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp \
       --gc_intervals ${TEST_DIR}/plasbin_utils_test_results/gc.txt \
       --out_file     ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd_gcp.csv \
       --log_file     ${TEST_DIR}/gc_probabilities_3.log

echo "----> Testing preprocessing with gene density and no GC probability"
python plasbin_utils.py preprocessing \
       --input_file   ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd.csv \
       --out_dir      ${TEST_DIR}/plasbin_utils_test_results_preprocessing \
       --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp_preprocessing \
       --gc_intervals ${TEST_DIR}/plasbin_utils_test_results/gc.txt \
       --out_file     ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd_preprocessing.csv \
       --log_file     ${TEST_DIR}/preprocessing_1.log

echo "----> Testing preprocessing with mappings and no gene density and no GC probability"
python plasbin_utils.py preprocessing \
       --input_file   ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c.csv \
       --out_dir      ${TEST_DIR}/plasbin_utils_test_results_preprocessing \
       --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp_preprocessing \
       --gc_intervals ${TEST_DIR}/plasbin_utils_test_results/gc.txt \
       --out_file     ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_preprocessing.csv \
       --log_file     ${TEST_DIR}/preprocessing_2.log

echo "----> Testing preprocessing with GC probability and no gene density and no mappings"
python plasbin_utils.py preprocessing \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_gcp.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results_preprocessing \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp_preprocessing \
       --gc_intervals  ${TEST_DIR}/plasbin_utils_test_results/gc.txt \
       --db_file       ${TEST_DIR}/plasbin_utils_test_results/pls.genes.fasta \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_gcp_preprocessing.csv \
       --log_file      ${TEST_DIR}/preprocessing_3.log

echo "----> Testing preprocessing with GC probability and mappings but no gene density"
python plasbin_utils.py preprocessing \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gcp.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results_preprocessing \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp_preprocessing \
       --gc_intervals  ${TEST_DIR}/plasbin_utils_test_results/gc.txt \
       --pid_threshold 0.95 \
       --cov_threshold 0.8 \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gcp_preprocessing.csv \
       --log_file      ${TEST_DIR}/preprocessing_4.log

echo "----> Testing preprocessing with GC probability and gene density and mappings (nothing done)"
python plasbin_utils.py preprocessing \
       --input_file    ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd_gcp.csv \
       --out_dir       ${TEST_DIR}/plasbin_utils_test_results_preprocessing \
       --tmp_dir       ${TEST_DIR}/plasbin_utils_test_tmp_preprocessing \
       --gc_intervals  ${TEST_DIR}/plasbin_utils_test_results/gc.txt \
       --out_file      ${TEST_DIR}/plasbin_utils_test_results/${INPUT_PREF}_g2c_gd_gcp_preprocessing.csv \
       --log_file      ${TEST_DIR}/preprocessing_5.log
