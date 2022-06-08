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

	contigs_dict[c]['Seed'] = 0		#Default
	contigs_dict[c]['PrPl'] = 0.01	#Default
	contigs_dict[c]['PrChr'] = 0.99	#Default
	contigs_dict[c]['log_ratio'] = math.log(0.01/0.99)	#Default

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
def get_data(assembly_file, contigs_dict, links_list):
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

#Reads the seed file and makes a list of seeds
def get_seeds(seeds_file, seeds_set):
	string_list = read_file(seeds_file)
	for line in string_list:
		line = line.split("\t")
		seeds_set.add(line[0])
	return seeds_set	

def get_class_probs(class_file, contigs_dict):
	string_list = read_file(class_file)
	for line in string_list[1:]:
		line = line.split(" ")
		c, prpl, prchr = line[0], float(line[3]), float(line[2])
		if prpl == 0:
			prpl = 0.01
		if prchr == 0:
			prchr = 0.01		
		contigs_dict[c]['PrPl'] = prpl
		contigs_dict[c]['PrChr'] = prchr
		contigs_dict[c]['log_ratio'] = math.log(prpl/prchr)

	return contigs_dict 

def get_gc_probs(gc_file, gc_probs):
	string_list = read_file(gc_file)
	for line in string_list:
		line = line.split("\t")
		c = line[0]
		gc_probs[c] = {}
		for i in range(1,len(line)):
			gc_probs[c][i] = float(line[i])
	return gc_probs

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
def contig_vars(m, contigs_dict, contigs, contigs_ext):
	contigs = {}
	contigs_ext = {}
	for c in contigs_dict:
		contigs[c] = m.addVar(vtype=GRB.BINARY, name='contig-'+c)
		contigs_ext[(c, 'h')] = m.addVar(vtype=GRB.BINARY, name='contigext-'+c+'h')
		contigs_ext[(c, 't')] = m.addVar(vtype=GRB.BINARY, name='contigext-'+c+'t')			
	return contigs, contigs_ext

#links[e] == 1 if link 'e' belongs to solution, else 0
#For each link e, we also add the edge variable in the reverse direction
def link_vars(m, links_list, contigs, links):
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
	bins = [1,2,3,4,5,6]
	for b in bins:
		plas_GC[b] = m.addVar(vtype=GRB.BINARY, name='plasmid-GC-bin-'+str(b))
	for c in gc_probs:
		contig_GC[c] = {}
		for b in bins:
			contig_GC[c][b] = m.addVar(vtype=GRB.BINARY, name='contig-'+c+'-GC-bin-'+str(b))
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

'''
#-------------------------------------------------------------
def rd_vars(m, contigs_dict, rd_diff, counted_rd_diff):
	for p in range(1):
		rd_diff[p] = {}
		counted_rd_diff[p] = {}
		for c in contigs_dict:	
			rd_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='rd-diff_contig-'+c+'_plasmid-'+str(p))
			counted_rd_diff[p][c] = m.addVar(vtype=GRB.CONTINUOUS, name='counted-rd-diff_contig-'+c+'_plasmid-'+str(p))
	return rd_diff, counted_rd_diff

	



#-----------------------------------------------------------------------		
'''





			


