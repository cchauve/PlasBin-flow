import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import gzip
from pathlib import Path
from Bio import SeqIO
import scipy.special as sc
import scipy.integrate as integrate
import numpy as np
from scipy.stats import binom
import math

def get_seq_details(seq, mol_type):
    '''
    Parsing entry in the fasta file and
    Updating the dictionary with all sequences
    '''
    sequence = str(seq.seq)
    gc_count = sequence.count('G') + sequence.count('C')
    length = len(sequence)
    return {
        'gc_count': gc_count,
        'length': length,
        'gc_percent': gc_count/length,
        'type': mol_type
    }

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

def compute_gc_intervals_files(
        chr_paths_file, pls_paths_file, out_txt_file, out_png_file
):
    '''
    Reads fasta sequences of plasmids and chromosomes in the reference database
    Computes and outputs details on gc content of molecules including
    violinplot of gc content distribution

    Args:
        - chr_paths_file (str): path to file containing the list of paths of chromosome FASTA files
        - pls_paths_file (str): path to file containing the list of paths of plasmids FASTA files
        - out_txt_file (str): written text file containing GC content of chromosomes and plasmids
        - out_png_file (str): written PNG file containing violin plot of GC content of chromosomes and plasmids
    
    Returns:
        None, creates files out_txt_file and out_png_file
    '''

    '''
    Dictionary with all sequences (plasmids as well as chromosomes)
    Key: Sequence Id, Values: Dictionary of attributes
    gc_count (int)
    length (int)
    gc_percent (float): Ratio of gc_count/length
    type (str): Type of molecule, plasmid or chromosome
    '''
    seq_dict = {}
    # Reading and storing input data for reference samples
    chr_paths = open(chr_paths_file,'r').read().splitlines()
    seq_dict = read_fastas(seq_dict, chr_paths, 'chromosome')
    pls_paths = open(pls_paths_file,'r').read().splitlines()
    seq_dict = read_fastas(seq_dict, pls_paths, 'plasmid')
    # Converting to dataframe
    seq_df = pd.DataFrame.from_dict(seq_dict).T
    seq_df.index.name='id'
    seq_df.to_csv(out_file)
    second_df = pd.DataFrame(seq_df.to_dict())
    # Violinplot from dataframe
    fig, ax = plt.subplots(figsize=(10,8))
    plt.suptitle('GC distribution')
    ax = sns.violinplot(x="type", y="gc_percent", data=second_df)
    plt.savefig(out_png_file)

''' Functions that should be extracted into a separate module '''

def read_file(filename):
    string = open(filename, "r").read()
    string_list = string.split("\n")
    string_list = [
        line
        for line in string_list
        if line and line[0] != '#'
    ]
    return string_list

# Storing contig details
## Stores the id of the contig
def get_id(line):
    return line[1]
## Stores the nucleotide sequence of the contig
def get_nucleotide_seq(line):
    return line[2]
## Stores the length of the sequence
def get_length(line):
    return int(line[3].split(':')[2])

# Takes a contig from the assembly file and initiates an entry in the contigs_dict
# Each contig is tagged with the following attributes:
# 1. Length of the contig (int)
# 2. Sequence of a contig (string)
def update_contigs_dict(contigs_dict, line):
    c = get_id(line)
    seq = get_nucleotide_seq(line)
    GC_cont = compute_GCratio(seq)
    n = get_length(line)
    
    contigs_dict[c] = {}
    contigs_dict[c]['Sequence'] = seq
    contigs_dict[c]['GC_cont'] = GC_cont
    contigs_dict[c]['Length'] = n
    
    return contigs_dict

# Reads the assembly file line by line and forwards a line
# to update_contigs_dict or get_link depending on the entry
def get_data(assembly_file, contigs_dict):
    string_list = read_file(assembly_file)
    for line in string_list:
        line = line.split("\t")
        if line[0] == 'S':
            contigs_dict = update_contigs_dict(contigs_dict, line)
    return contigs_dict

# Computes GC ratio: counts no. of 'G'/'C' occurences in the sequence and divide by the sequence length.
def compute_GCratio(seq):
    GC = 0
    ln_seq = 0
    for nucl in seq:
        if nucl == 'G' or nucl == 'C':
            GC += 1
        ln_seq += 1
    return GC/ln_seq

def combln(n, k):
    # Compute ln of n choose k. Note than n! = gamma(n+1).
    return sc.gammaln(n + 1) - sc.gammaln(k + 1) - sc.gammaln(n - k + 1)

def gprob2(n,g,p,m):
    # Compute probability of observing g GC nucleotides in a contig
    # of length n within a molecule of GC content p using pseducount m.
    # Done via logarithm to avoid overflow.
    alpha = m*p
    beta = m*(1-p)
    resultln = combln(n, g) + sc.betaln(g + alpha, n - g + beta) - sc.betaln(alpha, beta)
    return math.exp(resultln)

def compute_gc_probabilities_file(gfa_file, gc_intervals_file, out_file):

    if gc_intervals_file != None:
        probs = read_file(gc_intervals_file)
        probs = [float(x) for x in probs]
        probs = sorted(probs)
    else:
        probs = [0, 0.4,  0.45, 0.5, 0.55, 0.6, 1]

    gc_file = open(out_file, "w")
    gc_file.write('CTG')
    for i in range(0, len(probs)-1):
        gc_file.write('\t'+str(probs[i])+'-'+str(probs[i+1]))
    gc_file.write('\n')

    contigs_dict = {}
    contigs_dict = get_data(gfa_file, contigs_dict)
    
    m = 10
    for c in contigs_dict:
        n = contigs_dict[c]['Length']
        g_frac = contigs_dict[c]['GC_cont']
        g = g_frac*n
        total = 0
        gc_file.write(c)
        gp_array = []
        for i in range(0, len(probs)-1):
            gp2 = integrate.quad(lambda x: gprob2(n,g,x,m), probs[i], probs[i+1])
            gp2 = gp2[0]/(probs[i+1] - probs[i])
            total += gp2
            gp_array.append(gp2)
        for gp in gp_array:
            gc_file.write("\t"+str(gp / total))
        gc_file.write("\n")
