#!/usr/bin/python

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
	argparser = argparse.ArgumentParser()
	argparser.add_argument("contigs_fasta", help = "")
	argparser.add_argument("output_file", help = "")
	

	args = argparser.parse_args()	

	#Storing assembly contig sequences	
	contig_seqs = dict()
	contigs = set()
	with open(args.contigs_fasta) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			#print(record.id)
			contig_seqs[record.id] = record.seq
			contigs.add(record.id)
	
	component_seqs = dict()
	for contig in contig_seqs:
		component = contig.split("_")[7]
		if component not in component_seqs:
			component_seqs[component] = ''
		component_seqs[component] += contig_seqs[contig]
	
	record_list = []
	for comp in component_seqs:
		seq = component_seqs[comp]
		seq_len = len(seq)
		record = SeqRecord(Seq(seq),id='component_'+comp, description="length=%s"%(str(seq_len)))
		record_list.append(record)


    # write to file
	SeqIO.write(record_list, args.output_file, 'fasta')	
	#with open(args.output_file, "w") as output_handle:
	#	SeqIO.write(component_seqs, output_handle, "fasta")




main()