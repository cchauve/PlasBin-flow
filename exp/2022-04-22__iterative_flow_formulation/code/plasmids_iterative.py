__author__ = 'amane'

#-------------------

#USAGE: 
#time python plasmids_iterative.py --ag assembly.gfa --gc gc_probs.csv --c class_probs.csv --seeds seed_contigs.csv \
#				 --out output_dir --alpha1 alpha_1 --alpha2 alpha_2 --alpha3 alpha_3 --rmiter rmiter

from gurobipy import *
from sys import argv
import os
import time
from random import randint
import math
import argparse
import plasmids_preprocessing
import plasmids_postprocessing
import rmcircular

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

if __name__ == "__main__":
	#Parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--ag", help="Path to assembly graph file")
	parser.add_argument("--gc", help="Path to GC probabilities file")
	parser.add_argument("--c", help="Path to contig classification file")
	parser.add_argument("--seeds", help="Path to seed contigs file")
	parser.add_argument("--out", help="Path to output dir")
	parser.add_argument("--alpha1", nargs='?', const = 1, type=int, default = 1, help="Weight of flow term")
	parser.add_argument("--alpha2", nargs='?', const = 1, type=int, default = 1, help="Weight of GC content term")
	parser.add_argument("--alpha3", nargs='?', const = 1, type=int, default = 1, help="Weight of log probabilities term")	
	parser.add_argument("--rmiter", nargs='?', const = 1, type=int, default = 50, help="Number of iterations to remove circular components")

	args = parser.parse_args()

	output_dir = args.out
	assembly_file = args.ag
	gc_file = args.gc
	class_file = args.c
	seeds_file = args.seeds
	alpha1 = args.alpha1
	alpha2 = args.alpha2
	alpha3 = args.alpha3
	rmiter = args.rmiter

	#Naming and creating output files
	ratios = str(alpha1) + '.' + str(alpha2) + '.' + str(alpha3)
	output_folder = output_dir + '/' + ratios
	if not os.path.exists(output_folder):
	    os.makedirs(output_folder)	
		
	output_fasta = 'putative_plasmids.fasta'
	ques_fasta = 'questionable_plasmids.fasta'
	score_filename = 'MILP_objective.csv'
	output_contigs = 'contig_chains.csv'
	ques_contigs = 'questionable_contig_chains.csv'
	components = 'components.csv'

	#-----------------------------------------------
	#Main program
	contigs_dict = {}
	links_list = []
	seeds_set = set()
	gc_probs = {}
	capacities = {}

	contigs_dict, links_list = plasmids_preprocessing.get_data(assembly_file, contigs_dict, links_list)
	seeds_set = plasmids_preprocessing.get_seeds(seeds_file, seeds_set)
	contigs_dict = plasmids_preprocessing.get_class_probs(class_file, contigs_dict)
	gc_probs = plasmids_preprocessing.get_gc_probs(gc_file, gc_probs)
	#capacities = plasmids_preprocessing.get_caps(links_list, contigs_dict, capacities)

	components_dict = {}
	used_contigs = {}
	#ques_component_dict = {}
	
	len_total = 0
	n_contigs = 0
	for c in contigs_dict:
		if c in seeds_set:
			contigs_dict[c]['Seed'] = 1
		else:
			contigs_dict[c]['Seed'] = 0
		n_contigs += 1
		len_total += contigs_dict[c]['Length']		






	UBD_rd = 0
	UBD_GC = 0
	for c in contigs_dict:
		UBD_rd = max(UBD_rd, contigs_dict[c]['Read_depth'])
		UBD_GC = max(UBD_GC, contigs_dict[c]['GC_cont'])
	print(UBD_rd, UBD_GC, len_total)

	input_details = 'details.csv'
	details_file = open(os.path.join(output_folder, input_details), "w")
	details_file.write("Contig"+"\t"+"Read_depth"+ "\t"+"GC_cont"+"\t"+ "Length"+"\t"+ "PrPl"+"\t"+ "PrChr"+"\n")
	for c in contigs_dict:
		details_file.write(c+"\t"+str(contigs_dict[c]['Read_depth'])+ "\t"+str(contigs_dict[c]['GC_cont'])+"\t"+ str(contigs_dict[c]['Length'])+"\t"+ str(contigs_dict[c]['PrPl'])+"\t"+ str(contigs_dict[c]['PrChr'])+"\n")
		details_file.write("Link list:"+"\t")
		for link in links_list:
			if c == link[0][0] or c == link[1][0]:
				details_file.write(str(link)+"\t")
		details_file.write("\n")
				 
	output_fasta_file = open(os.path.join(output_folder, output_fasta), "w")
	ques_fasta_file = open(os.path.join(output_folder, ques_fasta), "w")
	score_file = open(os.path.join(output_folder, score_filename), "w")
	contigs_file = open(os.path.join(output_folder, output_contigs), "w")
	ques_contigs_file = open(os.path.join(output_folder, ques_contigs), "w")
	components_file = open(os.path.join(output_folder, components), "w")

	logfile_3 = open(os.path.join(output_folder,'circ_sequences.log'),"w")

	n_iter = 0
	q_iter = 0
	n_comp = 0
	q_comp = 0

	
	while len(seeds_set) > 0:
		#For consistency, a contig extremity should be a part of exactly one link for a specific plasmid.
		#Here, for each extremity, we make a list of links that involve the extremity.
		extr_dict = {}

		for c in contigs_dict:
			ext1, ext2 = (c, 'h'), (c, 't')
			extr_dict[ext1] = []
			extr_dict[ext2] = []

			for link in links_list:
				if link[0] == ext1 or link[1] == ext1:
					extr_dict[ext1].append(link) 
				if link[0] == ext2 or link[1] == ext2:
					extr_dict[ext2].append(link)

		incoming, outgoing = {}, {}
		for link in links_list:
			ext1, ext2 = link[0], link[1]

			if ext1 not in incoming:
				incoming[ext1] = [('S',ext1)]
				outgoing[ext1] = [(ext1,'T')]
			incoming[ext1].append(link[::-1])
			outgoing[ext1].append(link)

			if ext2 not in incoming:
				incoming[ext2] = [('S',ext2)]
				outgoing[ext2] = [(ext2,'T')]
			incoming[ext2].append(link)
			outgoing[ext2].append(link[::-1])		
		
		capacities = plasmids_preprocessing.get_caps(links_list, contigs_dict, capacities)
		print("\n\n\n\n\n")
		#-----------------------------------------------
		#Initializing the ILP
		m = Model("Plasmids")
		m.params.LogFile= os.path.join(output_folder,'m.log')
		m.setParam(GRB.Param.TimeLimit, 240.0)
		m.setParam(GRB.Param.MIPGap, 0.05)

		contigs = {}
		contigs_ext = {}
		contigs, contigs_ext = plasmids_preprocessing.contig_vars(m, contigs_dict, contigs, contigs_ext)

		links = {}
		links = plasmids_preprocessing.link_vars(m, links_list, links)

		plas_GC = {}
		contig_GC = {}
		plas_GC, contig_GC = plasmids_preprocessing.GC_vars(m, gc_probs, plas_GC, contig_GC)

		counted_seed = {}
		counted_seed = plasmids_preprocessing.seed_vars(m, contigs_dict, counted_seed)				

		#rd_diff, counted_rd_diff = {}, {}	
		#rd_diff, counted_rd_diff \
		#	= plasmids_preprocessing.rd_vars(m, contigs_dict, rd_diff, counted_rd_diff)

		counted_len = {}
		counted_len = plasmids_preprocessing.len_vars(m, contigs_dict, counted_len)

		flows = {}
		counted_flow = {}
		counted_overall_flow = {}
		flows, counted_flow, counted_overall_flow = plasmids_preprocessing.flow_vars(m, links, flows, counted_flow, counted_overall_flow)
		F = m.addVar(vtype=GRB.CONTINUOUS, name='overall-flow')	

		#print(links_list)
		print(len(links_list))
		#-----------------------------------------------
		#Setting up the expression for the objective function
		expr = LinExpr()
		alpha1, alpha2, alpha3 = float(alpha1), float(alpha2), float(alpha3)

		expr.addTerms(alpha1, F)
		for c in contigs:
			for b in plas_GC:
				expr.addTerms(-alpha2*(1-gc_probs[c][b]), contig_GC[c][b])
			expr.addTerms(alpha3*contigs_dict[c]['log_ratio'], contigs[c])
		m.setObjective(expr, GRB.MAXIMIZE)

		#-----------------------------------------------
		#Setting up constraints

		constraint_count = 0
		#Constraint type 1
		#A link 'e' is in 'p' only if both its endpoints are in 'p'

		for e in links:
			end1, end2 = e[0], e[1]
			if end1 != 'S' and end1 != 'T':
				m.addConstr(links[e] <= contigs_ext[end1], "link_ubd")
			if end2 != 'S' and end2 != 'T':	
				m.addConstr(links[e] <= contigs_ext[end2], "link_ubd")
			constraint_count += 2

		#Constraint type 2
		#Consistency: A contig extremity can occur only once in a plasmid
		#for extr in contigs_ext:
		#	c = extr[0]
		#	expr = LinExpr()
		#	for link in extr_dict[extr]:
		#		if link in links:
		#			expr.addTerms(1, links[link])
		#	m.addConstr(expr <= 1, "consistency")
		#	constraint_count += 1

		#Constraint type 2
		#An extremity is in 'p' only if at least one edge is incident on it
		for extr in contigs_ext:
			expr = LinExpr()
			link_count = 0
			for link in extr_dict[extr]:
				if link in links:
					link_count += 1
					expr.addTerms(1, links[link])
			m.addConstr(contigs_ext[extr] <= expr + 0.99, "extr_ubd")
			m.addConstr(contigs_ext[extr] >= expr, "extr_lbd")
			constraint_count += 2

		#Constraint type 3
		#A contig 'c' is in 'p' if at least one of its endpoints is in 'p'
		for c in contigs:
			end1, end2 = (c, 'h'), (c, 't')
			m.addConstr(contigs[c] >= contigs_ext[end1], "contig_lbd1")
			m.addConstr(contigs[c] >= contigs_ext[end2], "contig_lbd2")
			constraint_count += 2						
		
		#Constraint type 4
		#counted-length = length * existence-of-a-contig-in-plasmid
		for c in counted_len:
			if c in contigs:
				m.addConstr(counted_len[c] == contigs_dict[c]['Length']*contigs[c], "counted_length")
				constraint_count += 1
		expr = LinExpr()
		for c in contigs:
			expr.addTerms(1, counted_len[c])
		m.addConstr(expr <= 175000, "plasmid_length_upper_bound")
		constraint_count += 1	

		#Constraint type 5
		#counted-seed = seed * existence-of-a-contig-in-plasmid
		for c in counted_seed:
			if c in contigs:
				m.addConstr(counted_seed[c] == contigs_dict[c]['Seed']*contigs[c], "counted_seed")
				constraint_count += 1	

		#Constraint type 6
		#Computing flow out of S and choosing exactly one edge from S and exactly one into T 
		xS_expr = LinExpr()
		xT_expr = LinExpr()
		fS_expr = LinExpr()
		for e in links:
			if e[0] == 'S':
				xS_expr.addTerms(1, links[e])
				fS_expr.addTerms(1, counted_flow[e])
			if e[1] == 'T':
				xT_expr.addTerms(1, links[e])
		m.addConstr(xS_expr == 1, "one-edge-from-S")
		m.addConstr(xT_expr == 1, "one-edge-into-T")
		m.addConstr(F == fS_expr, "flow-from-S")

		#Constraint type 7
		#Flow conservation

		#Constraint type 8
		#Capacity constraints	
			
		for c in contigs:
			ext1, ext2 = (c, 'h'), (c, 't')
			hin_expr = LinExpr()
			hout_expr = LinExpr()
			tin_expr = LinExpr()
			tout_expr = LinExpr()		

			if ext1 in incoming:
				for e in incoming[ext1]:
					hin_expr.addTerms(1, counted_flow[e])
			if ext1 in outgoing:
				for e in outgoing[ext1]:
					hout_expr.addTerms(1, counted_flow[e])
			if ext2 in incoming:
				for e in incoming[ext2]:
					tin_expr.addTerms(1, counted_flow[e])
			if ext2 in outgoing:
				for e in outgoing[ext2]:
					tout_expr.addTerms(1, counted_flow[e])	

			m.addConstr(hin_expr == tout_expr, "contig-"+c+"-conservation-1")
			m.addConstr(tin_expr == hout_expr, "contig-"+c+"-conservation-2")

			m.addConstr(hin_expr + tin_expr <= contigs_dict[c]['Read_depth'], "cap-"+c)	

		for e in links:
			m.addConstr(flows[e] <= capacities[e], "cap-"+str(e))




		#Constraint type 9
		#Handling quadratic flow variables
		for e in links:
			m.addConstr(counted_flow[e] <= capacities[e]*links[e], "xf1")
			m.addConstr(counted_flow[e] <= flows[e], "xf2")
			m.addConstr(counted_flow[e] >= flows[e] - (1-links[e])*capacities[e], "xf3")
			m.addConstr(counted_flow[e] >= 0, "xf4")

			m.addConstr(counted_overall_flow[e] <= capacities[e]*links[e], "xF1")
			m.addConstr(counted_overall_flow[e] <= F, "xF2")
			m.addConstr(counted_overall_flow[e] >= F - (1-links[e])*capacities[e], "xF3")
			m.addConstr(counted_overall_flow[e] >= 0, "xF4")

		#Constraint type 10
		for c in contig_GC:
			for b in plas_GC:
				m.addConstr(contig_GC[c][b] <= contigs[c], "xGC1")
				m.addConstr(contig_GC[c][b] <= plas_GC[b], "xGC2")
				m.addConstr(contig_GC[c][b] >= contigs[c] + plas_GC[b] - 1, "xGC3")
				constraint_count += 3
		expr = LinExpr()
		for b in plas_GC:
			expr.addTerms(1, plas_GC[b])
		m.addConstr(expr == 1)		
					
		#Constraint type 11
		#Each plasmids should have at least one seed
		expr = LinExpr()
		for c in counted_seed:
			expr.addTerms(1, counted_seed[c])
		m.addConstr(expr >= 1, "seed_existence")
		constraint_count += 1

		#Constraint type 12
		#Preliminary path constraints
		c_expr = LinExpr()
		for c in contigs:
			c_expr.addTerms(1, contigs[c])
		l_expr = LinExpr()
		for e in links:
			l_expr.addTerms(1, links[e])
		m.addConstr(c_expr == l_expr - 1, "path_simple")

		#Constraint type 0
		#Non-negativity constraints
		for c in contigs:
			m.addConstr(contigs[c] >= 0, "non-negativity")
		for extr in contigs_ext:
			m.addConstr(contigs_ext[extr] >= 0, "non-negativity")
		for e in links:
			m.addConstr(links[e] >= 0, "non-negativity")

		for c in counted_len:
			m.addConstr(counted_len[c] >= 0, "non-negativity")
			m.addConstr(counted_seed[c] >= 0, "non-negativity")

		#TEST	
		start = time.time()
		m.optimize()
		stop = time.time()
		duration = stop - start

		if m.status == GRB.Status.INFEASIBLE:
			print ('The model cannot be solved because it is infeasible')
		elif m.status == GRB.Status.UNBOUNDED:
			print ('The model cannot be solved because it is unbounded')
		elif m.status == GRB.Status.INF_OR_UNBD:
			print ('The model cannot be solved because it is infeasible or unbounded ')

		if m.status == GRB.Status.INF_OR_UNBD or m.status == GRB.Status.INFEASIBLE or m.status == GRB.Status.UNBOUNDED:
			print ('The model cannot be solved because it is infeasible or unbounded ')
			m.computeIIS()
			m.write("m.ilp")
			for con in m.getConstrs():
				if con.IISConstr:
					print('%s' % con.constrName) 
			exit (1)	

		#Retrieving plasmid bin
		output_fasta_file = open(os.path.join(output_folder, output_fasta), "a")
		ques_fasta_file = open(os.path.join(output_folder, ques_fasta), "a")
		score_file = open(os.path.join(output_folder, score_filename), "a")
		contigs_file = open(os.path.join(output_folder, output_contigs), "a")
		ques_contigs_file = open(os.path.join(output_folder, ques_contigs), "a")
		components_file = open(os.path.join(output_folder, components), "a")

		#Flow zero condition
		if F.x == 0:
			exit (1)

		plasmid_length = 0

		for c in contigs:
			if contigs[c].x > 0:
				plasmid_length += contigs_dict[c]['Length']
		
		solution_links = set()
		solution_seq = ''
		soln_ext_dict = {}
		
		for e in links:
			if links[e].x > 0:
				end1, end2 = e[0], e[1]
				 
				#c1, ext1 = end1[0], end1[1]
				#c2, ext2 = end2[0], end2[1]			

				solution_links.add(e)
				if end1 not in soln_ext_dict:
					soln_ext_dict[end1] = [end2]
				else:
					soln_ext_dict[end1].append(end2)
				if end2 not in soln_ext_dict:
					soln_ext_dict[end2] = [end1]
				else:
					soln_ext_dict[end2].append(end1)	

		#Getting solution seq and contig chain			
		#solution_seq, contig_chain, plasmid_parts = plasmids_postprocessing.get_seq(solution_links, soln_ext_dict, contigs_dict, contigs)
		#Contig chain for single-contig plasmids
		#if len(contig_chain) == 0:
		#	for c in contigs:
		#		if contigs[c].x > 0:
		#			contig_chain.append(str(c)+'+')
						
		#contig_chain, plasmid_length, solution_seq = plasmids_postprocessing.refine_chain(contig_chain, contigs_dict)	

		#Sorting into putative and questionable
		if plasmid_length >= 1500:
			output_fasta_file.write(">plasmid_"+str(n_iter)+"\t"+"length="+str(plasmid_length)+"\n")
			#if len(solution_seq) != 0:							#if length of contig chain is more than zero
			#	output_fasta_file.write(solution_seq+"\n")
			#else:													#if plasmid is a single contig
			#	c = contig_chain[0][:-1]								
			#	output_fasta_file.write(contigs_dict[c]['Sequence']+"\n")
				
			#Printing putative contig chains to file
			#contigs_file.write("plasmid_"+str(n_iter)+";")
			#components_file.write("plasmid_"+str(n_iter)+";\n")
			#contigs_file.write(contig_chain[0])
			#if len(contig_chain) > 1:
			#	for x in contig_chain[1:]:
			#		contigs_file.write(","+x)					
			#contigs_file.write("\n")
			#if len(contig_chain) == 1:
			#	components_file.write(contig_chain[0])
			#	components_file.write("\n")
			#Printing putative contig chain details to file	
			contigs_file.write("#Contig\tPrPl\tFlow\tRead_depth\tLength\n")		
			#for x in contig_chain:
			#	c = x[:-1]
			for c in contigs:
				if contigs[c].x > 0.5:			
					contigs_file.write("# "+c+"\t"+str(contigs_dict[c]['PrPl'])+"\t"+str(F.x)+"\t"+str(contigs_dict[c]['Read_depth'])+"\t"+str(contigs_dict[c]['Length'])+"\n")				
			contigs_file.write("\n")

			contigs_file.write("#Link\tCapacity\n")
			for e in links:
				if links[e].x > 0:
					contigs_file.write("# "+str(e)+"\t"+str(capacities[e])+"\n")
			contigs_file.write("\n")

			#for x in plasmid_parts:
			#	components_dict[n_comp] = x
			#	n_comp += 1
			#	for i in range(len(x)):
			#		y = x[i]
			#		if i == 0:
			#			components_file.write(str(y[0]+y[1]))
			#		else:
			#			components_file.write(","+str(y[0]+y[1]))	
			#	components_file.write("\n")

			#Recording scores
			GC_sum = 0
			class_sum = 0
			for c in contigs:
				for b in contig_GC[c]:
					GC_sum += -alpha2*(1-gc_probs[c][b])*contig_GC[c][b].x
				class_sum += alpha3*contigs_dict[c]['log_ratio']*contigs[c].x
			score_file.write("putative_plasmid_"+str(n_iter)+"\t\t"+str(m.objVal)+"\t"+str(F.x)+"\t"+str(GC_sum)+"\t"+str(class_sum)+"\n")
			n_iter += 1					

		
		#Updating assembly graph and formulation
		for c in contigs:
			if contigs[c].x > 0:
				contigs_dict[c]['Read_depth'] = max(0, contigs_dict[c]['Read_depth'] - F.x)
				if c in seeds_set and contigs_dict[c]['Read_depth'] <= 0.5:
					contigs_dict[c]['Seed'] = 0
					seeds_set.remove(c)

				if contigs_dict[c]['Read_depth'] <= 0.05:
					#orig_c = c.split('_')[0]
					orig_c = c
					used_contigs[orig_c] = {}
					used_contigs[orig_c]['Length'] = contigs_dict[c]['Length']
					used_contigs[orig_c]['GC_cont'] = contigs_dict[c]['GC_cont']
					used_contigs[orig_c]['Sequence'] = contigs_dict[c]['Sequence']
					del(contigs_dict[c])	
					del(gc_probs[c])			

					temp_list = []
					for e in links_list:						
						if (c, 'h') not in e and (c, 't') not in e:
							temp_list.append(e)
					links_list = temp_list	
					
					
						
'''
		#Checking for existence of disconnected cyclic components	
		#Adding cycle elimination constraints if needed	
		logfile_3 = open(os.path.join(output_folder,'circ_sequences.log'),"a")
		circular = 1
		i = 1
		while circular == 1:
			start = time.time()
			m.optimize()
			stop = time.time()
			duration = stop - start

			if m.status == GRB.Status.INFEASIBLE:
				print ('The model cannot be solved because it is infeasible')
			elif m.status == GRB.Status.UNBOUNDED:
				print ('The model cannot be solved because it is unbounded')
			elif m.status == GRB.Status.INF_OR_UNBD:
				print ('The model cannot be solved because it is infeasible or unbounded ')

			if m.status == GRB.Status.INF_OR_UNBD or m.status == GRB.Status.INFEASIBLE or m.status == GRB.Status.UNBOUNDED:
				print ('The model cannot be solved because it is infeasible or unbounded ')
				m.computeIIS()
				m.write("m.ilp")
				for con in m.getConstrs():
					if con.IISConstr:
						print('%s' % con.constrName) 
				exit (1)

			solution_links = {}		
			solution_seq = {}
			soln_ext_dict = {}

			for p in links:
				solution_links[p] = set()
				solution_seq[p] = ''
				soln_ext_dict[p] = {}
				for e in links[p]:
					if links[p][e].x > 0:
						end1, end2 = e[0], e[1]
						#c1, ext1 = end1[0], end1[1]
						#c2, ext2 = end2[0], end2[1]			
						solution_links[p].add(e)
						if end1 not in soln_ext_dict[p]:
							soln_ext_dict[p][end1] = [end2]
						else:
							soln_ext_dict[p][end1].append(end2)
						if end2 not in soln_ext_dict[p]:
							soln_ext_dict[p][end2] = [end1]
						else:
							soln_ext_dict[p][end2].append(end1)	
				solution_seq[p], contig_chain, plasmid_parts = plasmids_postprocessing.get_seq(solution_links[p], soln_ext_dict[p], contigs_dict, contigs)			

			circ_seq = {}
			k = 1 
			logfile_3.write("Iteration "+ str(i)+"\n")

			circular = 0	

			circular, circ_seq, solution_seq[p] = \
				rmcircular.check_if_circular(solution_seq, circ_seq, k, i, soln_ext_dict, solution_links, contigs, logfile_3)
			
			if circular == 1:
				for x in circ_seq:
					chosen_seq = set()	
					chosen_seq = circ_seq[x]['Seq']
				
					expr = LinExpr()
					for e in chosen_seq:
						if e in links[p]:
							expr.addTerms(1, links[p][e])
						elif e[::-1] in links[p]:
							expr.addTerms(1, links[p][e[::-1]])
					logfile_3.write("\n\n")		
					m.addConstr(expr <= len(chosen_seq) - 1, "circular")
							
					i += 1
			if i >= rmiter:
				break
'''				
	#_______________
	#
	
'''	

			else:
				#Same for questionable plasmids		
				ques_fasta_file.write(">plasmid_"+str(n_iter)+"\t"+"length="+str(plasmid_length[p])+"\t"+"gene_density="+str(plasmid_gd[p])+"\t"+"mean_read_depth="+str(plasmid_rd[p])+"\n")
				if len(solution_seq[p]) != 0:
					ques_fasta_file.write(solution_seq[p]+"\n")
				else:
					c = contig_chain[0][:-1]
					ques_fasta_file.write(contigs_dict[c]['Sequence']+"\n")
		
				ques_contigs_file.write("plasmid_"+str(q_iter)+";")
				for x in contig_chain:
					if x == contig_chain[0]:
						ques_contigs_file.write(x)

					else:
						ques_contigs_file.write(","+x)
				ques_contigs_file.write("\n")		
				for x in contig_chain:
					c = x[:-1]
					ques_contigs_file.write("# "+c+"\t"+str(contigs_dict[c]['Gene_coverage'])+"\t"+str(contigs_dict[c]['Read_depth'])+"\t"+str(contigs_dict[c]['Length'])+"\n")				
				ques_contigs_file.write("\n")

				gd_sum, GC_sum = 0, 0
				for p in contigs:
					for c in contigs[p]:
						gd_sum += alpha1*(-contigs_dict[c]['Gene_coverage'])*contigs[p][c].x
						GC_sum += alpha2*counted_GC_diff[p][c].x
				score_file.write("questionable_plasmid_"+str(n_iter)+"\t"+str(m.objVal)+"\t"+str(gd_sum)+"\t"+str(GC_sum)+"\n")
				q_iter += 1
'''
'''	
	
		#Updating assembly graph and formulation
		for c in contigs:
			if contigs[c].x > 0:
				contigs_dict[c]['Read_depth'] = max(0, contigs_dict[c]['Read_depth'] - contigs[c].x)
				if c in seeds_set and contigs_dict[c]['Read_depth'] <= 0.5:
					contigs_dict[c]['Seed'] = 0
					seeds_set.remove(c)

				if contigs_dict[c]['Read_depth'] <= 0.05:
					orig_c = c.split('_')[0]
					used_contigs[orig_c] = {}
					used_contigs[orig_c]['Length'] = contigs_dict[c]['Length']
					used_contigs[orig_c]['GC_cont'] = contigs_dict[c]['GC_cont']
					used_contigs[orig_c]['Sequence'] = contigs_dict[c]['Sequence']
					del(contigs_dict[c])	
					del(gc_probs[c])			

					temp_list = []
					for e in links_list:						
						if (c, 'h') not in e and (c, 't') not in e:
							temp_list.append(e)
					links_list = temp_list
'''			
	#_______________
	#
'''
	#Retaining unique plasmids from output		
	
	component_contigs_list = []
	count = 0
	for k in components_dict: #For every individual plasmid component
		component = components_dict[k] 
		component_contigs = set([x[0].split('_')[0] for x in component])
		#(component_contigs)
		temp_list = []
		flag = 0 #flag = 1 if set of contigs of current plasmid is subset of set of contigs of any previous plasmids
		for contig_set in component_contigs_list:
			if contig_set[0].issubset(component_contigs) == False: 
				temp_list.append(contig_set)
			elif contig_set[0] == component_contigs:
				temp_list.append(contig_set)
			if flag == 0 and component_contigs.issubset(contig_set[0]):
				flag = 1
		if flag == 0:
			temp_list.append((component_contigs,count))
		component_contigs_list = temp_list	
		count += 1	
		
	for pair in component_contigs_list:
		idx = pair[1]
		temp_len = 0
		temp_seq = []
		gene_free_len = []
		gene_free_seq = []
		print(idx, component_contigs_list)
		chain = component_contigs_list[idx][0]
			
		for contig in chain:

			gd = used_contigs[contig]['Gene_coverage']
			ln = used_contigs[contig]['Length']

			if gd >= 0.3:
				if temp_len > 2000:
					gene_free_len.append(temp_len)
					gene_free_seq.append(temp_seq)
				temp_len = 0
				temp_seq = []
			else:
				temp_len += ln
				temp_seq.append(contig)

		for x in range(len(gene_free_seq)):
			components_list[idx] = plasmids_postprocessing.remove_sub_list(gene_free_seq,chain)			
			
	#Printing unique plasmids to output file
	output_fasta_file = open(os.path.join(output_folder, output_fasta), "w")
	ques_fasta_file = open(os.path.join(output_folder, ques_fasta), "w")
	components_file = open(os.path.join(output_folder, components), "w")
	count = 1
	for pair in component_contigs_list:
		idx = pair[1]
		chain = list(component_contigs_list[idx][0])
		if len(chain) > 1 and chain[0] == chain[-1]:
			chain = chain[:-1]
		length = 0
		gd = 0
		rd = 1
		seq = ''
		components_file.write("plasmid_"+str(count)+";")
		for contig in chain:
			length += used_contigs[contig]['Length']
			gd += used_contigs[contig]['Gene_coverage']*used_contigs[contig]['Length']
			seq += used_contigs[contig]['Sequence']
			components_file.write(str(contig)+",")
		gd = gd/length
		if gd >= 0.3 and length >= 1500: 	
			components_file.write("\n")

			output_fasta_file.write(">plasmid_"+str(count)+"\t"+"length="+str(length)+"\t"+"gene_density="+str(gd)+"\t"+"mean_read_depth="+str(1)+"\n")
			output_fasta_file.write(seq+"\n")
		else:
			ques_fasta_file.write(">plasmid_"+str(count)+"\t"+"length="+str(length)+"\t"+"gene_density="+str(gd)+"\t"+"mean_read_depth="+str(1)+"\n")
			ques_fasta_file.write(seq+"\n")			
		count+=1	
'''			
	
							


