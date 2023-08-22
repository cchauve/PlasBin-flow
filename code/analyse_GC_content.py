#!/usr/bin/python

#USAGE:  python analyse_GC_content.py --chr CHR_PATHS --pls PLS_PATHS --out OUT_FILE --vplot VPLOT

'''
Reads fasta sequences of plasmids and chromosomes in the reference database
Computes and outputs details on gc content of molecules including
Violinplot of gc content distribution
'''

import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import gzip
from pathlib import Path
from Bio import SeqIO

#Parsing arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("--chr", help = "File containing chromosome fasta paths")
argparser.add_argument("--pls", help = "File containing plasmid fasta paths")
argparser.add_argument("--out", help = "Path to output file")
argparser.add_argument("--vplot", help = "Path to violinplot output file")

args = argparser.parse_args()

#Paths to folders containing chromosome and plasmid fasta files (gunzipped)
CHR_PATHS = args.chr
PLS_PATHS = args.pls

#Output files
OUT_FILE = args.out
out_file = Path(OUT_FILE)
out_file.parent.mkdir(parents=True, exist_ok=True)

VPLOT = args.vplot
img_file = Path(VPLOT)
img_file.parent.mkdir(parents=True, exist_ok=True)

def get_seq_details(seq, mol_type):
	'''
	Parsing entry in the fasta file and
	Updating the dictionary with all sequences	
	'''
	sequence = str(seq.seq)
	gc_count = sequence.count('G') + sequence.count('C')
	length = len(sequence)
	return {'gc_count': gc_count, 'length': length, 'gc_percent': gc_count/length, 'type': mol_type}
	
def read_fastas(seq_dict, fasta_paths, mol_type):
	'''
	Read all paths to chromosome/plasmid fasta files
	Populate dictionary of sequences	
	'''
	for path in fasta_paths:
		with gzip.open(os.path.expanduser(path),'rt') as f:
			seqs = SeqIO.parse(f,'fasta')
			for seq in seqs:
				seq_id = seq.id
				seq_dict[seq_id] = get_seq_details(seq, mol_type)
	return seq_dict

'''
Dictionary with all sequences (plasmids as well as chromosomes)
Key: Sequence Id, Values: Dictionary of attributes
gc_count (int)
length (int)
gc_percent (float): Ratio of gc_count/length
type (str): Type of molecule, plasmid or chromosome 
'''
seq_dict = {}

#Reading and storing input data for reference samples
chr_paths = open(CHR_PATHS,'r').read().splitlines()
pls_paths = open(PLS_PATHS,'r').read().splitlines()

#Reading chromosome fasta files
seq_dict = read_fastas(seq_dict, chr_paths, 'chromosome')

#Reading plasmid fasta files		
seq_dict = read_fastas(seq_dict, pls_paths, 'plasmid')

#Converting to dataframe
seq_df = pd.DataFrame.from_dict(seq_dict).T
seq_df.index.name='id'
seq_df.to_csv(out_file)
second_df = pd.DataFrame(seq_df.to_dict()) #Adjusting dtypes from object to int/float as required

#Violinplot from dataframe
fig, ax = plt.subplots(figsize=(10,8))
plt.suptitle('GC distribution')
ax = sns.violinplot(x="type", y="gc_percent", data=second_df)
plt.savefig(img_file)