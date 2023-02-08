# PlasBin-flow
PlasBin-flow is a plasmid binning tool that uses a mixed integer linear programming model that relies on the concept of network flows. It accounts for several contig features inclduing sequencing coverage, the presence of plasmid genes and the GC content that often distinguishes plasmids from chromosomes to bin contigs into sets representing putative plasmids

## Overview
The code/ directory contains the source code of PlasBin-flow. The requirements to run PlasBin-flow have been listed below. The results/ directory contains the results of PlasBin, PlasBin-flow, plasmidSPAdes, MOB-recon, HyAsP and gplas.

## Requirements
The following are required to run PlasBin
1. Python (Version 3+; packages: networkx, random, math, sys)
2. Gurobi solver (Version 9.1.2+)

## Input
1. A file containing the details of the assembly graph (.gfa format), 
2. A file mapping genes from a plasmid marker database to contigs (.csv format)<br/>
Each line of this file contains the ids of the two sequences being mapped. In this case, we map genes to contigs, hence the file contains the gene name and the contig id. This is followed by the starting and ending positions of the contig sequence, denoting the region onto which the gene has been mapped.<br/>
The format for this file is the blastn output format 6. <br/>
3. A file containing a list of GC content likelihoods (one entry per line)<br/>
Each line of this file contains the probabilities that a particular contig belongs to a plasmid with GC content to various GC content ranges.

### Usage
```
python plasbin_flow.py --ag assembly.gfa --gc gc_probs.csv --map gene_contig_mapping.csv \
		--out output_dir --alpha1 alpha_1 --alpha2 alpha_2 --alpha3 alpha_3 --rmiter rmiter
```
Additional arguments
```
--rmiter			Number of iterations to remove circular components. (default: 50)
--alpha1			Weight of flow term. (default: 1)                              
--alpha2			Weight of GC content term. (default: 1)
--alpha3			Weight of gene density term. (default: 1)
```

## Output
### Components
List of components obtained in all iterations. </br>
File name: components.csv<br/>
Format: plasmid id;comma-separated list of contigs<br/>
Example: plasmid_0;23,25,10<br/>


### Unique components
The output of plasbin-flow may contain components the vertices of which form subsets of vertices in other components. In such cases, as a postprocessing step, we only consider the components with the largest vertex sets. We use the code unique_components.py to obtain such components. </br>
List of unique contig chains that are not subsets of other contig chains. </br>
File name: unique_components.csv<br/>
Format: plasmid id;comma-separated list of contigs<br/>
Example: plasmid_0;23,25,10<br/>
