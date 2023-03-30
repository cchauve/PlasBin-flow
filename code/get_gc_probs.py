import argparse
from sys import argv
import os
import scipy.special as sc
import scipy.integrate as integrate
#import scipy.stats
import numpy as np
from scipy.stats import binom
import math
import argparse

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

#Storing contig details
#-----------------------------------------------
#Stores the id of the contig
def get_id(line):
	return line[1]
#Stores the nucleotide sequence of the contig
def get_nucleotide_seq(line):
	#print(line[2])
	return line[2]		
#Computes GC ratio: counts no. of 'G'/'C' occurences in the sequence and divide by the sequence length.
def compute_GCratio(seq):
	GC = 0
	ln_seq = 0
	for nucl in seq:
		if nucl == 'G' or nucl == 'C':
			GC += 1
		ln_seq += 1
	return GC/ln_seq
#Stores the length of the sequence
def get_length(line):
	return int(line[3].split(':')[2])

#Takes a contig from the assembly file and initiates an entry in the contigs_dict
#Each contig is tagged with the following attributes:
#1. Length of the contig (int)
#2. Sequence of a contig (string)
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

#Reads the assembly file line by line and forwards a line 
#to update_contigs_dict or get_link depending on the entry
def get_data(assembly_file,contigs_dict):
	string_list = read_file(assembly_file)
	for line in string_list:
		line = line.split("\t")
		if line[0] == 'S':
			contigs_dict = update_contigs_dict(contigs_dict, line)
	return contigs_dict

def combln(n, k):
	#Compute ln of n choose k. Note than n! = gamma(n+1).
	return sc.gammaln(n + 1) - sc.gammaln(k + 1) - sc.gammaln(n - k + 1)

def gprob2(n,g,p,m):
	#Compute probability of observing g GC nucleotides in a contig
	#of length n within a molecule of GC content p using pseducount m. 
  	#Done via logarithm to avoid overflow.
	alpha = m*p
	beta = m*(1-p)
	resultln = combln(n, g) + sc.betaln(g + alpha, n - g + beta) - sc.betaln(alpha, beta)
	return math.exp(resultln)

if __name__ == "__main__":
	#Parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-ag", help="Path to assembly graph file")
	parser.add_argument("-outdir", help="Path to output dir")
	parser.add_argument("-outfile", help="Name of output file")
	parser.add_argument("-gcint", default = None, help="Path to GC interval file")

	args = parser.parse_args()

	output_dir = args.outdir
	output_file = args.outfile
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)


	assembly_file = args.ag
	gcint_file = args.gcint

	if gcint_file != None:
		probs = read_file(gcint_file)
		probs = [float(x) for x in probs]
		probs = sorted(probs)
	else:
		probs = [0, 0.4,  0.45, 0.5, 0.55, 0.6, 1]

	gc_file = open(os.path.join(output_dir, output_file), "w")
	gc_file.write('CTG')
	for i in range(0, len(probs)-1):
		gc_file.write('\t'+str(probs[i])+'-'+str(probs[i+1]))
	gc_file.write('\n')


	contigs_dict = {}
	contigs_dict = get_data(assembly_file, contigs_dict)
	
	
	m = 10 

	for c in contigs_dict:
		n = contigs_dict[c]['Length']
		g_frac = contigs_dict[c]['GC_cont']
		g = g_frac*n
		#print(g_frac)
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
		
