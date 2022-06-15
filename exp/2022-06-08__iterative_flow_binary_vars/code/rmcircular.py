from gurobipy import *
import collections
import random
import copy

def get_link_soln(soln_ext_dict, solution_links, links):
	for plasmid in links:	#Separating those with pva = 1
		soln_ext_dict[plasmid] = {}
		solution_links[plasmid] = []
		for link in links[plasmid]:
			var = links[plasmid][link]
			if var.x > 0:
				solution_links[plasmid].append(link) 
				soln_ext_dict[plasmid][link[0]] = link[1]
				soln_ext_dict[plasmid][link[1]] = link[0]
	return soln_ext_dict, solution_links

def other(ext):
	c = 't' if ext == 'h' else 'h'
	return c

def get_next_link(contig, ext, ext_dict):
	next_contig, next_ext, next_link = None, None, None
	other_ext = other(ext)[0]
	l = (contig, other_ext)
	if l in ext_dict:
		next_contig = ext_dict[l][0]
		next_ext = ext_dict[l][1]
	if next_contig != None:
		r = (next_contig, next_ext)
		next_link = (l,r)
	return next_link, next_contig, next_ext


def check_if_circular(solution_seq, circ_seq, k, i, neighbour_dict, soln_links, contigs, logfile_3):
	plasmid_links = copy.deepcopy(soln_links)
	n_dict = copy.deepcopy(neighbour_dict)

	extr_degree = {}
	plasmid_ends = {}
	#print(n_dict)
	for p in contigs:
		extr_degree[p] = {}
		plasmid_ends[p] = []
		solution_seq[p] = {}

		for c in contigs[p]:
			extr_degree[p][c] = {}
			if (c,'h') in n_dict[p]:
				extr_degree[p][c][(c,'h')] = len(n_dict[p][(c,'h')])
			else:
				extr_degree[p][c][(c,'h')] = 0
			if (c,'t') in n_dict[p]:		
				extr_degree[p][c][(c,'t')] = len(n_dict[p][(c,'t')])	
			else:
				extr_degree[p][c][(c,'t')] = 0				
	
	for p in extr_degree:
		for c in extr_degree[p]:
			if extr_degree[p][c][(c,'h')] > extr_degree[p][c][(c,'t')]:
				plasmid_ends[p].append((c,'h'))
			elif extr_degree[p][c][(c,'h')] < extr_degree[p][c][(c,'t')]:
				plasmid_ends[p].append((c,'t'))

	#print(extr_degree)
	#print(plasmid_ends)

	for p in plasmid_ends:
		if len(plasmid_ends[p]) == 0:
			linear_found = 1
		else:
			linear_found = 0	

		#print("LINEAR FOUND YET: ", linear_found)	

		seq_idx = 0	

		while len(plasmid_links[p]) != 0:
			plasmid_seg = collections.deque()	#Sequence of contig-orientation pairs
			current_seq = collections.deque()	#Sequence of links

			seq_idx += 1
			solution_seq[p][seq_idx] = {}

			if linear_found == 0:
				#Initial assignment
				#print("TRYING TO FIND LIN SEG")
				curr_vertex = plasmid_ends[p][0]
				c1 = curr_vertex[0]
				plasmid_seg.append((c1, '+')) if curr_vertex[1] == 'h' else plasmid_seg.append((c1, '-'))		
			
				while extr_degree[p][c1][curr_vertex] != 0: 
					#print(curr_vertex, n_dict[p][curr_vertex])

					#print(curr_vertex, n_dict)
					neighbour = n_dict[p][curr_vertex][0]
					n_dict[p][curr_vertex].remove(neighbour)
					n_dict[p][neighbour].remove(curr_vertex)
					#print(curr_vertex, neighbour, n_dict[curr_vertex])
					if len(n_dict[p][curr_vertex]) == 0:
						del(n_dict[p][curr_vertex])
					if len(n_dict[p][neighbour]) == 0:
						del(n_dict[p][neighbour])	

					curr_link = (curr_vertex, neighbour)
					if curr_link not in plasmid_links[p]:
						curr_link = curr_link[::-1]
				
					#print(curr_link)
					#print(plasmid_links)
					plasmid_links[p].remove(curr_link)


				
					ext1, ext2 = curr_link[0], curr_link[1]
					extr_degree[p][ext1[0]][ext1] -= 1
					extr_degree[p][ext2[0]][ext2] -= 1

					c1 = neighbour[0]
					plasmid_seg.append((c1, '+')) if neighbour[1] == 't' else plasmid_seg.append((c1, '-'))
					curr_vertex = (c1, other(neighbour[1]))

				#print(plasmid_seg)		
				solution_seq[p][seq_idx]['Seq'] = plasmid_seg
				solution_seq[p][seq_idx]['Type'] = 'L'
				linear_found = 1	


			else:
				#print("TRYING TO FIND CIRCULAR SEG")
				circ_seq[k] = {}
				#print(len(n_dict), n_dict)
				curr_vertex = random.choice(list(n_dict[p]))		
				c1 = curr_vertex[0]
				plasmid_seg.append((c1, '+')) if curr_vertex[1] == 'h' else plasmid_seg.append((c1, '-'))
				#plasmid_seg.append((c2, '+')) if neighbour[1] == 't' else plasmid_seg.append((c2, '-'))			
				
				while extr_degree[p][c1][curr_vertex] != 0: 
					#print(curr_vertex, n_dict[p][curr_vertex])

					#print(curr_vertex, n_dict)
					neighbour = n_dict[p][curr_vertex][0]
					n_dict[p][curr_vertex].remove(neighbour)
					n_dict[p][neighbour].remove(curr_vertex)
					#print(curr_vertex, neighbour, n_dict[curr_vertex])
					if len(n_dict[p][curr_vertex]) == 0:
						del(n_dict[p][curr_vertex])
					if len(n_dict[p][neighbour]) == 0:
						del(n_dict[p][neighbour])	

					curr_link = (curr_vertex, neighbour)
					if curr_link not in plasmid_links[p]:
						curr_link = curr_link[::-1]
					
					current_seq.append(curr_link)
					#print(curr_link)
					#print(plasmid_links)
					plasmid_links[p].remove(curr_link)


					
					ext1, ext2 = curr_link[0], curr_link[1]
					extr_degree[p][ext1[0]][ext1] -= 1
					extr_degree[p][ext2[0]][ext2] -= 1

					c1 = neighbour[0]
					plasmid_seg.append((c1, '+')) if neighbour[1] == 't' else plasmid_seg.append((c1, '-'))
					curr_vertex = (c1, other(neighbour[1]))

				circ_seq[k]['Plasmid'] = p
				circ_seq[k]['Seq'] = current_seq	
				logfile_3.write(str(i)+"\t"+str(current_seq)+"\n")
				k += 1
				#print(current_seq)	
				#plasmid.append(plasmid_seg)	
				solution_seq[p][seq_idx]['Seq'] = plasmid_seg
				solution_seq[p][seq_idx]['Type'] = 'C'	

	#print("INSIDE RMCIRCULAR ", circ_seq)						

	if len(circ_seq) == 0:
		circular = 0
	else:
		circular = 1				
	return circular, circ_seq, solution_seq			