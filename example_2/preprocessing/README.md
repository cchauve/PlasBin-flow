# Preprocessing samples for PlasBin-flow

The file `input.csv` contains an example (composed of three
*Campylobacter* samples) of running PlasBin-flow preprocessing.

To preprocess the three samples, using the GC intervals computed in
the tuning (see [../tuning/](../tuning)) the command is
```
python ../../code/plasbin_utils.py preprocessing \
       --input_file    input.csv \
       --out_dir       output \
       --tmp_dir       preprocessing_tmp \
       --log_file      preprocessing.log \
       --out_file      preprocessing_output.csv \
       --db_file       ../tuning/output/pls.genes.fasta \
       --pid_threshold 0.9 \
       --cov_threshold 0.7 \
       --gc_intervals  ../tuning/output/gc.txt
```

The optional parameters `pid_threshold 0.9` is non-default (default
value = 0.95), as are `--cov_threshold 0.7` (default value = 0.8) and
`--gc_intervals ../tuning/output/gc.txt` (default value as in file
[../default/gc_intervals.txt](../default/gc_intervals.txt)).

The file `preprocessing_output.csv` will contain fields for the
contigs plasmid score and GC content probabilities files for all three
samples.

The temporary directory `preprocessing_tmp` is deleted but if the
optional parameter `--keep_tmp_dir` is also used.

The main results will be the files
- `output/<sample>.gc.tsv`: GC content probabilities file;
- `output/<sample>.gd.tsv`: plasmid scores (gene density) file.

Additionally, for each sample, the following files will be created:
- `<sample>.genes_mappings.tsv`: mapping of reference plasmid
  genes to sample contigs;
- `<sample>.gd.tsv`: gene density file (used as contigs plasmid score
  file);
- `<sample>.ground_truth.tsv`: ground truth file.

The file `preprocessing.log` contains a detailed log.
