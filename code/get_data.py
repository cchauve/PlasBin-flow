from __future__ import division
import math

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
#Stores the read depth of the contig
def get_read_depth(line):
	return float(line[4].split(':')[2])		

#Takes a contig from the assembly file and initiates an entry in the contigs_dict
#Each contig is tagged with the following attributes:
#1. Length of the contig (int)
#2. Overall read depth of the contig (float)
#3. GC content of the contig (float)
#4. Indication if the contig is a seed (binary)
#5. Sequence of a contig (string)
#6. Probability that a contig is of plasmidic origin (float)
#7. Probability that a contig is of chromosomal origin (float)
#8. Log ratio of the above probabilities (float)
def update_contigs_dict(contigs_dict, line):
	c = get_id(line)
	seq = get_nucleotide_seq(line) 
	GC_cont = compute_GCratio(seq)
	n = get_length(line)
	rd = get_read_depth(line)

	contigs_dict[c] = {}
	contigs_dict[c]['Sequence'] = seq
	contigs_dict[c]['GC_cont'] = GC_cont
	contigs_dict[c]['Length'] = n
	contigs_dict[c]['Read_depth'] = rd
	contigs_dict[c]['Gene_coverage'] = 0				#Default
	contigs_dict[c]['Gene_coverage_intervals'] = []		#Default

	return contigs_dict

#A link is of the type: ((l1, e1),(l2, e2)) 
#where l1, l2 are adjacent links and e1, e2 are the connected link extremities
def get_link(line, links_list):
	c1, o1, c2, o2 = line[1], line[2], line[3], line[4]
	if o1 == '+':
		ext1 = 'h'
	else:
		ext1 = 't'
	if o2 == '+':
		ext2 = 't'
	else:
		ext2 = 'h'

	edge = ((c1, ext1),(c2, ext2))
	links_list.append(edge)
	
	return links_list

#Reads the assembly file line by line and forwards a line 
#to update_contigs_dict or get_link depending on the entry
def get_ag_details(assembly_file,contigs_dict, links_list):
	string_list = read_file(assembly_file)
	count_s = 0
	count_l = 0
	for line in string_list:
		line = line.split("\t")
		if line[0] == 'S':
			contigs_dict = update_contigs_dict(contigs_dict, line)
		elif line[0] == 'L':
			links_list = get_link(line, links_list)		
	return contigs_dict, links_list

#Determine seeds from class probs
def get_seeds(contigs_dict, seeds_set):
	for c in contigs_dict:
		if contigs_dict[c]['Length'] >= 2650:
			#if (contigs_dict[c]['Gene_coverage'] >= 0.8 or contigs_dict[c]['PrPl'] >= 0.8) \
			#	  and contigs_dict[c]['PrChr'] <= 0.2:
			if contigs_dict[c]['Gene_coverage'] >= 0.58:
				seeds_set.add(c)
	return seeds_set

#Takes the gene covering intervals for a contig and finds their union
#The length of the union is used to compute gene coverage
def get_union(intervals):
	union = []
	for begin, end in sorted(intervals):
		if union and union[-1][1] >= begin - 1:
			union[-1][1] = max(union[-1][1], end)
		else:
			union.append([begin, end])
	return union		

#Computes the gene coverage for each contig
def get_gene_coverage(mapping_file, contigs_dict):
	string_list = read_file(mapping_file)
	possible_seeds = []
	for line in string_list:
		line = line.split("\t")	
		qseqid, sseqid = line[0], line[1]
		sstart, send = line[8], line[9]

		if sseqid in contigs_dict:
			if sseqid not in possible_seeds:
				possible_seeds.append(sseqid)
			if int(sstart) > int(send):
				contigs_dict[sseqid]['Gene_coverage_intervals'].append((int(send), int(sstart)))
			else:
				contigs_dict[sseqid]['Gene_coverage_intervals'].append((int(sstart), int(send)))

	for sseqid in contigs_dict:
		union = get_union(contigs_dict[sseqid]['Gene_coverage_intervals'])
		length = contigs_dict[sseqid]['Length']
		covered = 0
		for interval in union:
			covered += interval[1] - interval[0] + 1
		contigs_dict[sseqid]['Gene_coverage'] = covered/length
	return contigs_dict	

def get_gc_probs(gc_file, gc_probs, gc_pens):
	string_list = read_file(gc_file)
	gc_intervals = string_list[0].split("\t")[1:]
	for line in string_list[1:]:
		line = line.split("\t")
		c = line[0]
		gc_probs[c] = {}
		gc_pens[c] = {}
		max_prob = 0
		for i in range(len(gc_intervals)):
			gc_int = gc_intervals[i]
			bin_prob = float(line[i+1])
			gc_probs[c][gc_int] = bin_prob 
			if max_prob <= bin_prob:
				max_prob = bin_prob
		for gc_int in gc_probs[c]:
			gc_pens[c][gc_int] = gc_probs[c][gc_int] - max_prob
	return gc_probs, gc_pens

def get_caps(links_list, contigs_dict, capacities):
	for e in links_list:
		c1, c2 = e[0][0], e[1][0]
		cap = min(contigs_dict[c1]['Read_depth'], contigs_dict[c2]['Read_depth'])
		capacities[e] = cap
		capacities[e[::-1]] = cap

	for c in contigs_dict:
		ext1, ext2 = (c, 'h'), (c, 't')
		cap = contigs_dict[c]['Read_depth']
		capacities[('S', ext1)] = cap
		capacities[('S', ext2)] = cap
		capacities[(ext1, 'T')] = cap
		capacities[(ext2, 'T')] = cap		
	return capacities