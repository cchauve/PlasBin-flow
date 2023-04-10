#!/usr/bin/python

# Analyses how well predicted plasmids resp. their contigs cover reference plasmids.
# It determines the proportions of reference plasmids covered by individual / all predicted plasmids (and vice versa)
# and uses these information to score the predictions (recall, precision, F1 score).
from __future__ import division
import argparse
import pandas as pd
from Bio import SeqIO

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

# read BLAST output file into table
def read_blast_output(file):
	col_names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]  # outfmt 6
	return pd.read_csv(file, sep = '\t', names = col_names, dtype = str)

# compute number of positions (integers) covered by a list of potentially overlapping intervals
def num_covered_positions(intervals):
	intervals.sort(key = lambda x: x[0])  # intervals is now sorted by start position

	num_pos_covered = 0
	last_pos_covered = 0  # last (right-most) position of contig covered so far
	for start, end in intervals:
		if end <= last_pos_covered:
			pass  # contained in previous interval -> no new position covered
		else:
			num_pos_covered += end - max(last_pos_covered + 1, start) + 1
			last_pos_covered = end

	return num_pos_covered

def analyse_mapping(blast_out, contig_seqs, plasmid_seqs, contigs, plasmids, out):
	covered_pred_sections = dict([(ctg, []) for ctg in contig_seqs])
	covered_pred_sections_per_ref = dict()

	for ctg in contig_seqs:
		for pls in sorted(plasmid_seqs):
			ctg_pls_hits = blast_out.loc[blast_out.qseqid == ctg].loc[blast_out.sseqid == pls]
			covered_pred_sections_per_ref[(ctg, pls)] = []

			for index, row in ctg_pls_hits.iterrows():
				qstart = int(row[6])
				qend = int(row[7])
				interval = (qstart, qend) if qstart <= qend else (qend, qstart)
				covered_pred_sections[ctg].append(interval)
				covered_pred_sections_per_ref[(ctg, pls)].append(interval)

	#Write to TSV file
	out.write("#PLS\tCTG\tMAP\tPLS_LEN\tCTG_LEN\n")
	for pls in sorted(plasmid_seqs):
		for ctg in sorted(contig_seqs):
			val = num_covered_positions(covered_pred_sections_per_ref[(ctg, pls)]) / len(contig_seqs[ctg])
			ctg_len = len(contig_seqs[ctg])
			pls_len = len(plasmid_seqs[pls])
			out.write(pls + '\t' + ctg + '\t' + str(val) + '\t' + str(pls_len) + '\t' + str(ctg_len) + '\n') 
		out.write("\n")

	#out.write(ref_pls+'\t'+bin_order[count]+'\t'+str(val)+'\t'+str(lengths[ref_pls])+'\t'+str(lengths[bin_order[count]])+'\n')

	



def main():
	argparser = argparse.ArgumentParser()
	argparser.add_argument("plasmids_fasta", help = "")
	argparser.add_argument("contigs_fasta", help = "")
	argparser.add_argument("blast_results", help = "")
	argparser.add_argument("output_file", help = "")
	

	args = argparser.parse_args()	
    
	#Storing true / predicted plasmid sequences
	plasmid_seqs = dict()
	plasmids = set()
	with open(args.plasmids_fasta) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			#print(record.id)
			plasmid_seqs[record.id] = record.seq
			plasmids.add(record.id)

	#Storing assembly contig sequences	
	contig_seqs = dict()
	contigs = set()
	with open(args.contigs_fasta) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			#print(record.id)
			contig_seqs[record.id] = record.seq
			contigs.add(record.id)

	blast_out = read_blast_output(args.blast_results)

	with open(args.output_file, "w") as out:
		analyse_mapping(blast_out, contig_seqs, plasmid_seqs, contigs, plasmids, out)
	


main()

