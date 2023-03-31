# PlasBin-flow usage example

Here we show how to use PlasBin-flow to extract sets of contigs from an assembly graph to form plasmid bins. We provide the commands to generate the input for PlasBin-flow and executing PlasBin-flow itself.

##  Creating a gene database from a collection of plasmids

As a first step, we create the database of genes that is to be used for mapping onto contigs. Note that if you already have a database of genes in FASTA format, this step can be skipped. This step requires the use of standard UNIX tools (curl and rm).
```
python ../code/get_gd.py create ../database/genes.fasta -a ../database/accessions.txt
```
The list of accession ids used for building the gene database used in our experiments has been provided in the file `database/accessions.txt`. In this example, we store the gene sequences in `database/genes.fasta`.

## Mapping genes to contigs

We then map the genes onto contigs. Note that this step requires the use of UNIX command rm as well as BLAST+ suite of tools (v2.6.0; makeblastdb and blastn).
```
python ../code/get_gd.py map ../database/genes.fasta input/gene_contig_mapping.tsv -ag input/assembly.gfa
```
For this example, we use the sequences of plasmids genes in the database, stored in the file `database/genes.fasta` generated in the previous step. The assembly graph file required for this step is `assembly.gfa` under the `input` folder. We store the output (the mapping of genes to contigs in blastn output format 6)in the file `gene_contig_mapping.tsv` in the same folder. 

## Creating the GC content file

PlasBin-flow requires, for each contig, the probablities that the contig originates from a molecule of GC content within some pre-defined GC content ranges. We run the following command to compute these probabilities.
```
python ../code/get_gc_probs.py -ag input/assembly.gfa -outdir input -outfile gc_probs.tsv -gcint `input/gc_intervals.txt`
```
Here, `input/assembly.gfa` is the assembly graph file. The list of endpoints of preferred GC content ranges can also be provided as an optional input, with one endpoint per line. This list should include `0` and `1` as mandatory endpoints. All other endpoints should be numbers between `0` and `1`. The file `input/gc_intervals.txt` lists the endpoints that are used to determine the GC content ranges in this example.

The generated output is stored in `input/gc_probs.tsv`. The first column of the file lists the contig ids while the remaining list the probabilities that a contig originates from a molecule of GC content within the respective GC content range defined by using the endpoints.

## Usage

All the files required as input for PlasBin-flow are provided in the `input` folder. We now use the fpllowing command to execute PlasBin-flow.
```
python ../code/plasbin_flow.py -ag input/assembly.gfa -gc input/gc_probs.tsv \
		-map input/gene_contig_mapping.tsv -outdir output -outfile test.out
```
where `input/assembly.gfa` is the assembly graph file, `input/gc_probs.tsv` is the GC content file described above and `input/gene_contig_mapping.tsv` is the gene mapping file described above. The output directory specified is `output` while the output file `test.out` contains the plasmid bins identified by PlasBin-flow. The output format has been described below. 

Additional arguments:
```
-rmiter			Maximum number of iterations to remove circular components. (optional, default: 50)
-alpha1			Weight of flow term. (optional, default: 1)                              
-alpha2			Weight of GC content term. (optional, default: 1)
-alpha3			Weight of gene density term. (optional, default: 1)
```

## Output
The output of PlasBin-flow is a TSV file (`output/test.out` for this example) with each line containing the following information:
```
Plasmid bin			Number (ID) associated with the plasmid bin
Flow value			Flow value associated with the plasmid bin
GC bin				Index of GC content interval associated with the plasmid bin
Contigs				Comma-separated list of contigs associated with plasmid bin along with their multiplicities
```





