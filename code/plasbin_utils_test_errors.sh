#!/bin/bash
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --account=def-chauvec
#SBATCH --output=/home/chauvec/projects/ctb-chauvec/PLASMIDS/USRA_PLASMIDS_2023/bin/PlasBin-flow/code/plasbin_utils_test_errors.out
#SBATCH --error=/home/chauvec/projects/ctb-chauvec/PLASMIDS/USRA_PLASMIDS_2023/bin/PlasBin-flow/code/plasbin_utils_test_errors.out
#SBATCH --job-name=plasbin_utils_test_errors

module load StdEnv/2020  gcc/9.3.0 blast+/2.12.0
module load python/3

TEST_DIR=/home/chauvec/projects/ctb-chauvec/PLASMIDS/USRA_PLASMIDS_2023/bin/PlasBin-flow/code/plasbin_utils_test

echo "--> Test with incorrect data"

rm -rf ${TEST_DIR}/plasbin_utils_test_results_bug* ${TEST_DIR}/plasbin_utils_test_tmp_bug*

echo "----> Dataset file with no header"
python plasbin_utils.py pls_genes_db \
       --input_file ${TEST_DIR}/plasbin_utils_test_bug1.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug1 \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug1 \
       --log_file   ${TEST_DIR}/plasbin_utils_test_bug1.log
echo "----> Dataset file with incorrect header"
python plasbin_utils.py pls_genes_db \
       --input_file ${TEST_DIR}/plasbin_utils_test_bug2.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug2 \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug2 \
       --log_file   ${TEST_DIR}/plasbin_utils_test_bug2.log
echo "----> Missing file entry in dataset file"
python plasbin_utils.py pls_genes_db \
       --input_file ${TEST_DIR}/plasbin_utils_test_bug3.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug3 \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug3 \
       --log_file   ${TEST_DIR}/plasbin_utils_test_bug3.log
echo "----> Missing ground truth field in dataset file"
python plasbin_utils.py seeds \
       --input_file ${TEST_DIR}/plasbin_utils_test.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug4a \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug4a \
       --db_file    ${TEST_DIR}/plasbin_utils_test_results/pls.genes.fasta \
       --log_file   ${TEST_DIR}/plasbin_utils_test_bug4a.log
python plasbin_utils.py tuning \
       --input_file ${TEST_DIR}/plasbin_utils_test.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug4b \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug4b \
       --log_file   ${TEST_DIR}/plasbin_utils_test_bug4b.log
echo "-------> Non-existing file"
python plasbin_utils.py pls_genes_db \
       --input_file ${TEST_DIR}/plasbin_utils_test_bug4.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug7a \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug7a \
       --log_file   ${TEST_DIR}/plasbin_utils_test_bug4c.log
python plasbin_utils.py gc_intervals \
       --input_file ${TEST_DIR}/plasbin_utils_test_bug4.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug7b \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug7b \
       --log_file   ${TEST_DIR}/plasbin_utils_test_bug4d.log
echo "----> Missing plasmid genes database"
python plasbin_utils.py map_genes_to_ctgs \
       --input_file ${TEST_DIR}/plasbin_utils_test.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug5a \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug5a \
       --db_file    ${TEST_DIR}/plasbin_utils_test_results_bug5a/xxpls.genes.fasta \
       --log_file   ${TEST_DIR}/plasbin_utils_test_bug5a.log
python plasbin_utils.py seeds \
       --input_file ${TEST_DIR}/plasbin_utils_test_gt.csv \
       --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug5b \
       --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug5b \
       --db_file    ${TEST_DIR}/plasbin_utils_test_results5b/xxpls.genes.fasta \
       --log_file   ${TEST_DIR}/plasbin_utils_test_bug5b.log
echo "----> Missing GC intervals file"
python plasbin_utils.py gc_probabilities \
       --input_file   ${TEST_DIR}/plasbin_utils_test.csv \
       --out_dir      ${TEST_DIR}/plasbin_utils_test_results_bug6e \
       --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp_bug6e \
       --gc_intervals ${TEST_DIR}/xxgc_intervals.txt \
       --out_file     ${TEST_DIR}/plasbin_utils_test_gcp5c.csv \
       --log_file     ${TEST_DIR}/plasbin_utils_test_bug5c.log

# echo "----> Corrupted dataset/sequence file"

# echo "----> Incorrect GC intervals file"
# python plasbin_utils.py gc_probabilities \
#        --input_file   ${TEST_DIR}/plasbin_utils_test.csv \
#        --out_dir      ${TEST_DIR}/plasbin_utils_test_results_bug6a \
#        --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp_bug6a \
#        --gc_intervals ${TEST_DIR}/gc_intervals_bug1.txt \
#        --out_file     ${TEST_DIR}/plasbin_utils_test_gcp6a.csv \
#        --log_file     ${TEST_DIR}/plasbin_utils_test_bug6a.log
# python plasbin_utils.py gc_probabilities \
#        --input_file   ${TEST_DIR}/plasbin_utils_test.csv \
#        --out_dir      ${TEST_DIR}/plasbin_utils_test_results_bug6b \
#        --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp_bug6b \
#        --gc_intervals ${TEST_DIR}/gc_intervals_bug2.txt \
#        --out_file     ${TEST_DIR}/plasbin_utils_test_gcp6b.csv \
#        --log_file     ${TEST_DIR}/plasbin_utils_test_bug6b.log
# python plasbin_utils.py gc_probabilities \
#        --input_file   ${TEST_DIR}/plasbin_utils_test.csv \
#        --out_dir      ${TEST_DIR}/plasbin_utils_test_results_bug6c \
#        --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp_bug6c \
#        --gc_intervals ${TEST_DIR}/gc_intervals_bug3.txt \
#        --out_file     ${TEST_DIR}/plasbin_utils_test_gcp6c.csv \
#        --log_file     ${TEST_DIR}/plasbin_utils_test_bug6c.log
# python plasbin_utils.py gc_probabilities \
#        --input_file   ${TEST_DIR}/plasbin_utils_test.csv \
#        --out_dir      ${TEST_DIR}/plasbin_utils_test_results_bug6d \
#        --tmp_dir      ${TEST_DIR}/plasbin_utils_test_tmp_bug6d \
#        --gc_intervals ${TEST_DIR}/gc_intervals_bug4.txt \
#        --out_file     ${TEST_DIR}/plasbin_utils_test_gcp6d.csv \
#        --log_file     ${TEST_DIR}/plasbin_utils_test_bug6d.log
# TO DO: Incorrect ground truth

# echo "-------> Incorrect entry in gfa file"
# python plasbin_utils.py map_genes_to_ctgs \
#        --input_file ${TEST_DIR}/plasbin_utils_test_bug5.csv \
#        --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug8 \
#        --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug8 \
#        --db_file    ${TEST_DIR}/plasbin_utils_test_results/pls.genes.fasta \
#        --log_file   ${TEST_DIR}/plasbin_utils_test_bug8.log
# echo "-------> Incorrect entry(ctg/plasmid/chromosome) in pls/pls/chr file"
# python plasbin_utils.py pls_genes_db \
#        --input_file ${TEST_DIR}/plasbin_utils_test_bug6a.csv \
#        --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug9a \
#        --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug9a \
#        --log_file   ${TEST_DIR}/plasbin_utils_test_bug9a.log
# python plasbin_utils.py gc_intervals \
#        --input_file ${TEST_DIR}/plasbin_utils_test_bug6a.csv \
#        --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug9b \
#        --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug9b \
#        --log_file   ${TEST_DIR}/plasbin_utils_test_bug9b.log
# python plasbin_utils.py gc_intervals \
#        --input_file ${TEST_DIR}/plasbin_utils_test_bug6b.csv \
#        --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug9c \
#        --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug9c \
#        --log_file   ${TEST_DIR}/plasbin_utils_test_bug9c.log
# python plasbin_utils.py gc_intervals \
#        --input_file ${TEST_DIR}/plasbin_utils_test_bug6c.csv \
#        --out_dir    ${TEST_DIR}/plasbin_utils_test_results_bug9d \
#        --tmp_dir    ${TEST_DIR}/plasbin_utils_test_tmp_bug9dx \
#        --log_file   ${TEST_DIR}/plasbin_utils_test_bug9d.log

# # TO DO: non-DNA sequences in data files
