#!/bin/bash
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --account=AACOUNT_NAME
#SBATCH --output=plasbin-flow_%A_%a.out
#SBATCH --error=plasbin-flow_%A_%a.err
#SBATCH --job-name=plasbin-flow
#SBATCH --array=1-2

source ${PLASBIN_FLOW_VIRTUAL_ENV}/bin/activate

INPUT=input.csv
PREPROCESSING_DIR=../preprocessing
GC_INTERVALS=../tuning/output/gc.txt

sed '1d' ${PREPROCESSING_DIR}/preprocessing_output.csv > ${INPUT}

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${INPUT} | cut -f1 -d',')
ASSEMBLER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${INPUT} | cut -f2 -d ',')
GFA=${PREPROCESSING_DIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${INPUT} | cut -f3 -d ',')
GC=${PREPROCESSING_DIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${INPUT} | cut -f4 -d ',')
GD=${PREPROCESSING_DIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${INPUT} | cut -f6 -d ',')

echo ${SLURM_ARRAY_TASK_ID} ${SAMPLE} 

rm -rf output/${SAMPLE}*

python ../../code/plasbin_flow.py \
       -ag ${GFA} \
       -gc ${GC} \
       -score ${GD} \
       -out_dir output \
       -out_file ${SAMPLE}.pred.txt \
       -log_file ${SAMPLE}.log \
       -assembler ${ASSEMBLER} \
       -gc_intervals ${GC_INTERVALS}
done
