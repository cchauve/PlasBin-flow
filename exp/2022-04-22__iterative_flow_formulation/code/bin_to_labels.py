from sys import argv
import os
import argparse

import plasmids_preprocessing

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line] #Read line only if it is nonempty
	return string_list

if __name__ == "__main__":
	#Parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--ag", help="Path to assembly graph file")
	parser.add_argument("--bin", help="Path to predicted bins file")
	parser.add_argument("--out", help="Path to output file")
	
	args = parser.parse_args()

	output_file = args.out
	assembly_file = args.ag
	bin_file = args.bin

	contigs_dict = {}
	links_list = []
	#seeds_set = set()
	#gc_probs = {}
	#capacities = {}

	contigs_dict, links_list = plasmids_preprocessing.get_data(assembly_file, contigs_dict, links_list)
	
	string_list = read_file(bin_file)
	bin_dict = {}

	delimit = 0
	for line in string_list:
		line = line.split("\t")
		#print(line[0])
		if "Contig" in line[0] or "Link" in line[0]:
			delimit += 1
		else:	 
			if delimit < 2:
				c = line[0][2:]
				bin_dict[c] = 2

	for c in contigs_dict:
		if c not in bin_dict:
			bin_dict[c] = 1			

	labels_file = open(output_file, "w")
	labels_file.write("short_read_contig_id,short_read_contig_length,bin_id,hybrid_contig_length,number_of_residue_matches,is_circular,label\n")

	for c in contigs_dict:
		sr_id = c
		sr_len = contigs_dict[c]['Length']
		bin_id = bin_dict[c]
		bin_len = 0
		res_mat = 0
		circ = 0
		if bin_id == 1:
			label = 'chromosome'
		else:
			label = 'plasmid'	
		labels_file.write(c+","+str(sr_len)+","+str(bin_id)+","+str(bin_len)+","+str(res_mat)+",False,"+label+"\n")
