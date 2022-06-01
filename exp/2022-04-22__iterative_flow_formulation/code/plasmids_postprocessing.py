from gurobipy import *
import collections
import random
import copy

def other(ext):
	c = 't' if ext == 'h' else 'h'
	return c	

def get_seq(soln_links, neighbour_dict, contigs_dict, contigs):
	plasmid_links = copy.deepcopy(soln_links)
	n_dict = copy.deepcopy(neighbour_dict)

	reached_links = {}
	plasmid_parts = []

	extr_degree = {}
	for c in contigs:
		extr_degree[c] = {}
		if (c,'h') in n_dict:
			extr_degree[c][(c,'h')] = len(n_dict[(c,'h')])
		else:
			extr_degree[c][(c,'h')] = 0
		if (c,'t') in n_dict:		
			extr_degree[c][(c,'t')] = len(n_dict[(c,'t')])	
		else:
			extr_degree[c][(c,'t')] = 0
	extr_degree['S'] = {}
	extr_degree['S']['S'] = len(n_dict['S'])	

	extr_degree['T'] = {}
	extr_degree['T']['T'] = len(n_dict['T'])		


			

	plasmid_ends = []
	for c in extr_degree:
		if c == 'S' or c == 'T':
			plasmid_ends.append(c)
		else:	
			if extr_degree[c][(c,'h')] > extr_degree[c][(c,'t')]:
				plasmid_ends.append((c,'h'))
			elif extr_degree[c][(c,'h')] < extr_degree[c][(c,'t')]:
				plasmid_ends.append((c,'t'))

	if len(plasmid_ends) == 0:
		linear_found = 1
	else:
		linear_found = 0	

	while len(plasmid_links) != 0:
		plasmid_seg = collections.deque()
		if linear_found == 1:
			#Initial assignment
			curr_vertex = plasmid_ends[0]
			c1 = curr_vertex[0]
			if c1 == 'S' or c1 == 'T':
				plasmid_seg.append(c1)
			else:	
				plasmid_seg.append((c1, '+')) if curr_vertex[1] == 'h' else plasmid_seg.append((c1, '-'))			
			
			while extr_degree[c1][curr_vertex] != 0: 
				neighbour = n_dict[curr_vertex][0]
				n_dict[curr_vertex].remove(neighbour)
				n_dict[neighbour].remove(curr_vertex)

				if len(n_dict[curr_vertex]) == 0:
					del(n_dict[curr_vertex])
				if len(n_dict[neighbour]) == 0:
					del(n_dict[neighbour])	

				curr_link = (curr_vertex, neighbour)
				if curr_link not in plasmid_links:
					curr_link = curr_link[::-1]
				
				plasmid_links.remove(curr_link)

				ext1, ext2 = curr_link[0], curr_link[1]
				if ext1 != 'S' and ext1 != 'T':
					extr_degree[ext1[0]][ext1] -= 1
				if ext2 != 'S' and ext2 != 'T':	
					extr_degree[ext2[0]][ext2] -= 1

				c1 = neighbour[0]
				if c1 == 'S' or c1 == 'T':
					plasmid_seg.append(c1)
				else:	
					plasmid_seg.append((c1, '+')) if neighbour[1] == 't' else plasmid_seg.append((c1, '-'))

			plasmid_parts.append(plasmid_seg)	
			linear_found = 1	

		else:
			curr_vertex = random.choice(list(n_dict))		
			c1 = curr_vertex[0]
			if c1 == 'S' or c1 == 'T':
				plasmid_seg.append(c1)
			else:	
				plasmid_seg.append((c1, '+')) if curr_vertex[1] == 'h' else plasmid_seg.append((c1, '-'))
			
			while extr_degree[c1][curr_vertex] != 0: 
				neighbour = n_dict[curr_vertex][0]
				n_dict[curr_vertex].remove(neighbour)
				n_dict[neighbour].remove(curr_vertex)
				if len(n_dict[curr_vertex]) == 0:
					del(n_dict[curr_vertex])
				if len(n_dict[neighbour]) == 0:
					del(n_dict[neighbour])	

				curr_link = (curr_vertex, neighbour)
				if curr_link not in plasmid_links:
					curr_link = curr_link[::-1]
				
				plasmid_links.remove(curr_link)

				ext1, ext2 = curr_link[0], curr_link[1]
				if ext1 != 'S' and ext1 != 'T':
					extr_degree[ext1[0]][ext1] -= 1
				if ext2 != 'S' and ext2 != 'T':	
					extr_degree[ext2[0]][ext2] -= 1

				c1 = neighbour[0]
				if c1 == 'S' or c1 == 'T':
					plasmid_seg.append(c1)
				else:	
					plasmid_seg.append((c1, '+')) if neighbour[1] == 't' else plasmid_seg.append((c1, '-'))
				
				print(c1, neighbour)
				if c1 == 'S' or c1 == 'T':
					curr_vertex = (neighbour)
				else:	
					curr_vertex = (c1, other(neighbour[1]))

			plasmid_parts.append(plasmid_seg)
	plasmid_seq = ''
	contig_chain = []						
	for seg in plasmid_parts:
		for contig in seg:
			c, sign = contig[0], contig[1]
			if c in contigs_dict:
				if sign == '+':
					contig_chain.append(str(c)+'+')
					plasmid_seq += contigs_dict[c]['Sequence']
				else:
					plasmid_seq += contigs_dict[c]['Sequence'][::-1]	
					contig_chain.append(str(c)+'-')				
	return plasmid_seq, contig_chain, plasmid_parts
		
def trim_zero_end(contig_chain, contigs_dict, end_pt):
	l , r = 0, 0
	if end_pt == 'l':
		#print(contig_chain)
		c = contig_chain[0][:-1]
		prpl = contigs_dict[c]['PrPl']
		if prpl == 0.01:
			contig_chain = contig_chain[1:]
			l = 0
		else:
			l = 1
	else:
		c = contig_chain[-1][:-1]
		prpl = contigs_dict[c]['PrPl']
		if prpl == 0:
			contig_chain = contig_chain[:-1]
			r = 0
		else:
			r = 1						
	return(contig_chain, l, r)

def get_rev_compliment(string):
	rev_string = string.replace("A","-").replace("T","A").replace("-","T")
	rev_string = rev_string.replace("C","-").replace("G","C").replace("-","G")
	return rev_string			

def refine_chain(contig_chain, contigs_dict):
	l, r = 0, 0
	end_pt = 'l'
	while l == 0:
		contig_chain, l, r = trim_zero_end(contig_chain, contigs_dict, end_pt)
	end_pt = 'r'
	while r == 0:
		contig_chain, l, r = trim_zero_end(contig_chain, contigs_dict, end_pt)	

	length, rd = 0, 0
	seq = ''
	for x in contig_chain:
		c, sign = x[:-1], x[-1]
		length += contigs_dict[c]['Length']
		if contigs_dict[c]['Read_depth'] < 1:
			rd += contigs_dict[c]['Read_depth']*contigs_dict[c]['Length']
		else:
			rd += 1*contigs_dict[c]['Length']
		if sign == '+':	
			seq += contigs_dict[c]['Sequence']	
		else:
			rev_seq = get_rev_compliment(contigs_dict[c]['Sequence'][::-1])
			seq += rev_seq
	rd = rd/length				
	return contig_chain, length, seq

# Returns that starting and ending point (index) of the sublist, if it exists, otherwise 'None'.
def find_sub_list(sub_list, in_list):
	sub_list_length = len(sub_list[0])
	flag = 0
	for i in range(len(in_list)-sub_list_length):
		if sub_list[0] == in_list[i:i+sub_list_length]:
			flag = 1
			return (i,i+sub_list_length)
	if flag == 0:
		return None	       

# Removes the sublist, if it exists and returns a new list, otherwise returns the old list.
def remove_sub_list(sub_list, in_list):
	indices = find_sub_list(sub_list, in_list)
	if not indices is None:
		return in_list[0:indices[0]] + in_list[indices[1]:]
	else:
		return in_list
		