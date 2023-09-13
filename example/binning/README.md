# Preprocessing and processing samples for PlasBin-flow

The file `input.csv` contains an example (composed of two
assemblies of a *Pseudomonas aeruginosa* sample, one Unicycler
assembly and one SKESA assembly) of running PlasBin-flow
preprocessing and binning commands.

To preprocess the two samples, using the default GC intervals and the
default plasmids gene database, the command is
```
python ../../code/plasbin_utils.py preprocessing \
       --input_file    input.csv \
       --out_dir       preprocessing_output \
       --tmp_dir       preprocessing_tmp \
       --log_file      preprocessing.log \
       --out_file      preprocessing_output.csv \
       --db_file       ../../database/genes.fasta \
       --pid_threshold 0.9 \
       --cov_threshold 0.7
```

The optional parameters `pid_threshold 0.9` is non-default (default
value = 0.95), as is `--cov_threshold 0.7` (default value = 0.8.

The file `preprocessing_output.csv` will contain fields for the
contigs plasmid score and GC content probabilities files for both
samples.

The temporary directory `preprocessing_tmp` is deleted but if the
optional parameter `--keep_tmp_dir` is also used.

The main results will be the files
- `preprocessing_output/<sample>.gc.tsv`: GC content probabilities file;
- `preprocessing_output/<sample>.gd.tsv`: plasmid scores (gene density) file.

Additionally, for each sample, the following files will be created:
- `<sample>.genes_mappings.tsv`: mapping of reference plasmid
  genes to sample contigs;
- `<sample>.gd.tsv`: gene density file (used as contigs plasmid score
  file);
- `<sample>.ground_truth.tsv`: ground truth file.

The file `preprocessing.log` contains a detailed log.

After the preprocessing is completed, the samples can be processed to
compute plasmid bins, using the default GC intervals, using the command
```
for TASK_ID in {1..2};
do
    SAMPLE_LINE=`expr ${TASK_ID} + 1`
    SAMPLE=$(sed -n "${SAMPLE_LINE}p" preprocessing_output.csv | cut -f1 -d',')
    ASSEMBLER=$(sed -n "${SAMPLE_LINE}p" preprocessing_output.csv | cut -f2 -d ',')
    GFA=$(sed -n "${SAMPLE_LINE}p" preprocessing_output.csv | cut -f3 -d ',')
    GC=$(sed -n "${SAMPLE_LINE}p" preprocessing_output.csv | cut -f4 -d ',')
    GD=$(sed -n "${SAMPLE_LINE}p" preprocessing_output.csv | cut -f6 -d ',')

    python ../../code/plasbin_flow.py \
           -ag ${GFA} \
           -gc ${GC} \
           -score ${GD} \
           -out_dir binning_output \
           -out_file ${SAMPLE}.pred.txt \
           -log_file binning_${SAMPLE}.log \
           -assembler ${ASSEMBLER}
done
```
This command uses a BASH loop and process both samples iteratively,
but samples could be processed on a HPC system, using for example the
array feature of the SLURM scheduler.

The output of the binning command are the files
```
binning_output/SAMEA12292472-skesa.pred.txt
binning_output/SAMEA12292472-unicycler.pred.txt
```
that each contains the predicted bins for the input
assembly graph.