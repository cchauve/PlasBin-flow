# Tuning PlasBin-flow

The file `input.csv` contains a toy example (composed of five
*Campylobacter* samples, for proper tuning a large number of samples
is required) of the data required to tune PlasBin-flow parameters.

To tune PlasBin-flow using this tuning dataset an the reference
plasmid genes database provided with PlasBin-flow and non-default
parameters, we first create a reference plasmid genes database:

```
python ../../code/plasbin_utils.py pls_genes_db \
       --input_file input.csv \
       --out_dir output \
       --tmp_dir tuning_tmp
```
which creates a file `output/pls.genes.fasta`.

Then we run the tuning command:

```
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
```

The optional parameters `pid_threshold 0.9` is non-default (default
value = 0.95), as are `--cov_threshold 0.7` (default value = 0.8) and
`--n_gcints n_gcints 5` (default value = 6).

The optional parameters `--out_file tuning_output.csv` and `--log_file
tuning.log` will result in creating an updated dataset CSV file and a
detailed log file.

The temporary directory `tuning_tmp` is deleted but if the optional
parameter `--keep_tmp_dir` is also used.

The main results will be the files
- `output/gc.[csv,png,txt]` where `output/gc.txt` is the GC intervals
  file;
- `output/seeds.txt` that contains the optimal seed contigs thresholds
  for the considered tuning dataset.

Note that the GC intervals file `output/gc.txt` shows very small
intervals due to the fact the tuning was done on a very small set of
very similar samples, where the GC content of chromosomes and plasmids
is quite similar, which can be seen in the files `output/gc.csv` ad
`output/gc.png`.  The seeds parameters file `output/seeds.txt` shows
many optimal seed parameters due to the small size of the reference
plasmid genes database.

Additionally, for each sample, the following files will be created:
- `<sample>.genes_mappings.tsv`: mapping of reference plasmid
  genes to sample contigs;
- `<sample>.gd.tsv`: gene density file (used as contigs plasmid score
  file);
- `<sample>.ground_truth.tsv`: ground truth file.