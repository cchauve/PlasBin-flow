from __future__ import division
from gurobipy import *
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
	contigs_dict[c]['Seed'] = 0		#Default
	contigs_dict[c]['PrPl'] = 0.5	#Default
	contigs_dict[c]['PrChr'] = 0.5	#Default
	contigs_dict[c]['log_ratio'] = math.log(0.5/0.5)	#Default

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
def get_data(assembly_file,contigs_dict, links_list):
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

def get_class_probs(class_file, contigs_dict):
	string_list = read_file(class_file)
	for line in string_list[1:]:
		#print(line)
		prpl, prchr = 0, 0

		#PlasForest results
		line = line.split(",")
		c = line[1]
		prpl = float(line[2])
		prchr = float(line[3])

		'''
		#mlplasmids results
		line = line.split(" ")
		c, prpl, prchr = line[0], float(line[3]), float(line[2])
		if prpl == 0:
			prpl = 0.01
		if prchr == 0:
			prchr = 0.01
		'''			
		if c in contigs_dict:
			contigs_dict[c]['PrPl'] = prpl
			contigs_dict[c]['PrChr'] = prchr
			#contigs_dict[c]['log_ratio'] = math.log(prpl/prchr)

	return contigs_dict 

def get_log_ratio(contigs_dict, p):
	for c in contigs_dict:
		prpl = contigs_dict[c]['PrPl']
		prchr = contigs_dict[c]['PrChr']
		gd = contigs_dict[c]['Gene_coverage']

		contigs_dict[c]['log_ratio'] = math.log((p*prpl+(1-p)*gd)/prchr)
	return contigs_dict

def get_gc_probs(gc_file, gc_probs, gc_pens):
	string_list = read_file(gc_file)
	for line in string_list:
		line = line.split("\t")
		c = line[0]
		gc_probs[c] = {}
		gc_pens[c] = {}
		max_prob = 0
		for i in range(1,len(line)):
			gc_probs[c][i] = float(line[i])
			if max_prob <= float(line[i]):
				max_prob = float(line[i])
		for i in range(1,len(line)):
			gc_pens[c][i] = gc_probs[c][i] - max_prob
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

#-----------------------------------------------------------	
#contigs[c] == 1 if contig 'c' belongs to solution, else 0
def contig_vars(m, contigs_dict, contigs):
	contigs = {}
	for c in contigs_dict:
		contigs[c] = m.addVar(vtype=GRB.BINARY, name='contig-'+c)		
	return contigs

#links[e] == 1 if link 'e' belongs to solution, else 0
#For each link e, we also add the edge variable in the reverse direction
def link_vars(m, links_list, links, contigs):
	links = {}
	for e in links_list:		
		links[e] = m.addVar(vtype=GRB.BINARY, name='link-'+e[0][0]+e[0][1]+e[1][0]+e[1][1])
		links[e[::-1]] = m.addVar(vtype=GRB.BINARY, name='link-'+e[1][0]+e[1][1]+e[0][0]+e[0][1])

	for c in contigs:
		h_ext = (c, 'h')
		t_ext = (c, 't')
		
		links[('S',h_ext)] = m.addVar(vtype=GRB.BINARY, name='link-S-to-'+c+'-h')
		links[(h_ext,'T')] = m.addVar(vtype=GRB.BINARY, name='link-'+c+'-h-to-T')
		links[('S',t_ext)] = m.addVar(vtype=GRB.BINARY, name='link-S-to-'+c+'-t')
		links[(t_ext,'T')] = m.addVar(vtype=GRB.BINARY, name='link-'+c+'-t-to-T')			
	return links

#TODO
def GC_vars(m, gc_probs, plas_GC, contig_GC):
	#bins = [1,2,3,4,5,6]

	#for b in bins:
	#	plas_GC[b] = m.addVar(vtype=GRB.BINARY, name='plasmid-GC-bin-'+str(b))
	for c in gc_probs:
		contig_GC[c] = {}
		for b in gc_probs[c]:
			contig_GC[c][b] = m.addVar(vtype=GRB.BINARY, name='contig-'+c+'-GC-bin-'+str(b))
			if b not in plas_GC:
				plas_GC[b] = m.addVar(vtype=GRB.BINARY, name='plasmid-GC-bin-'+str(b))
	return plas_GC, contig_GC		

#TODO
def seed_vars(m, contigs_dict, counted_seed):
	for c in contigs_dict:		
		counted_seed[c] = m.addVar(vtype=GRB.BINARY, name='counted-seed_contig-'+c)
	return counted_seed

#TODO
def len_vars(m, contigs_dict, counted_len):
	counted_ln = {}	
	for c in contigs_dict:
		counted_ln[c] = m.addVar(vtype=GRB.INTEGER, name='counted-ln_contig-'+c)
	return counted_ln

#TODO
def flow_vars(m, links, flows, counted_overall_flow):
	for e in links:
		flows[e] = m.addVar(vtype=GRB.CONTINUOUS, name='flow_edge-'+str(e[0])+'-to-'+str(e[1]))	
		counted_overall_flow[e] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-overall-flow_edge-'+str(e[0])+'-to-'+str(e[1]))	
	return flows, counted_overall_flow	







			


