__author__ = 'amane'

#-------------------

#USAGE: 
#time python plasmids_flow.py --ag assembly.gfa --gc gc_probs.csv --c class_probs.csv --seeds seed_contigs.csv \
#				 --out output_dir --alpha1 alpha_1 --alpha2 alpha_2 --alpha3 alpha_3 --rmiter rmiter

from re import L
from gurobipy import *
from sys import argv
import os
import time
from random import randint
import math
import csv
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
	var_vals = 'var_values.txt'

	#-----------------------------------------------
	#Main program
	contigs_dict = {}   #Key: contig IDs, Values: Contig attributes (provided as or derived from input)
	links_list = []     #List of unordered links: Each link is a pair of extrmeities (e.g. (('1','h'),('2','t')))
	seeds_set = set()   #Set of contig IDs from seeds file
	gc_probs = {}       #Key: contigs IDs, Values: Probability for each GC bin
	capacities = {}     #Key: 

	contigs_dict, links_list = plasmids_preprocessing.get_data(assembly_file, contigs_dict, links_list)
	seeds_set = plasmids_preprocessing.get_seeds(seeds_file, seeds_set)
	contigs_dict = plasmids_preprocessing.get_class_probs(class_file, contigs_dict)
	gc_probs = plasmids_preprocessing.get_gc_probs(gc_file, gc_probs)

	for c in contigs_dict:  #Assigning seeds
		if c in seeds_set:
			contigs_dict[c]['Seed'] = 1
		else:
			contigs_dict[c]['Seed'] = 0

	#Keeping a log of details about contigs in the assembly graph
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
	var_file = open(os.path.join(output_folder, var_vals), "w")

	logfile_3 = open(os.path.join(output_folder,'circ_sequences.log'),"w")	

	n_iter = 0
	q_iter = 0
	n_comp = 0
	q_comp = 0

	while len(seeds_set) > 0:
		#For consistency, both extremities for a contig should be part of exactly the same number of links in a plasmid.
		#For each extremity, we make a list of links that involve the extremity.
		extr_dict = {}	#Key: Extremity (e.g. ('1','h')), Value: List of unordered links incident on extremity
		#List of incoming and outgoing edges for each extremity
		#Including links from S and links to T
		incoming, outgoing = {}, {}	#Key: Extremity, Value: List of ordered links in/out of extremity

		for c in contigs_dict:
			ext1, ext2 = (c, 'h'), (c, 't')
			extr_dict[ext1] = [('S',(c, 'h')),((c, 'h'),'T')]
			extr_dict[ext2] = [('S',(c, 't')),((c, 't'),'T')]

			incoming[ext1] = [('S', ext1)]
			outgoing[ext1] = [(ext1, 'T')]
			incoming[ext2] = [('S', ext2)]
			outgoing[ext2] = [(ext2, 'T')]

			for link in links_list:
				if link[0] == ext1 or link[1] == ext1:
					extr_dict[ext1].append(link) 
					extr_dict[ext1].append(link[::-1]) 
				if link[0] == ext2 or link[1] == ext2:
					extr_dict[ext2].append(link)
					extr_dict[ext2].append(link[::-1])

		#print("EXTR", extr_dict)

		UBD_rd, LBD_rd = 0, 100
		for c in contigs_dict:	#Computing the upper bound for the read depth (and hence the upper bound for the flow)
			UBD_rd = max(UBD_rd, contigs_dict[c]['Read_depth'])
			LBD_rd = min(LBD_rd, contigs_dict[c]['Read_depth'])
		print(UBD_rd, LBD_rd)		

		#TODO: Include in above loop.	
		for link in links_list:
			ext1, ext2 = link[0], link[1]
			incoming[ext1].append(link[::-1])
			outgoing[ext1].append(link)
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

		#Initializing variables
		contigs = {}	#Key: Contig (e.g. '1'), Value: Gurobi binary variable
		#ext = {}		#Key: Extremity (e.g. ('1','h')), Value: Gurobi binary variable
		contigs = plasmids_preprocessing.contig_vars(m, contigs_dict, contigs)

		links = {}		#Key: Directed link from one extremity to another 
						#(e.g. (('1','h'),('2','t')) ), Value: Gurobi binary variable
		links = plasmids_preprocessing.link_vars(m, links_list, contigs, links)

		plas_GC = {}	#Key: GC bin, Value: Gurobi binary variable
		contig_GC = {}	#Nested Dict: Key: Contig, Value: Dict of GC bins (Key: GC bin, Value: Gurobi binary variable)
		plas_GC, contig_GC = plasmids_preprocessing.GC_vars(m, gc_probs, plas_GC, contig_GC)

		#counted_seed = {}	#Key: Contig, Value: Gurobi binary variable
		#counted_seed = plasmids_preprocessing.seed_vars(m, contigs_dict, counted_seed)				

		#counted_len = {}	#Key: Contig, Value: Gurobi integer variable
		#counted_len = plasmids_preprocessing.len_vars(m, contigs_dict, counted_len)

		flows = {}		#Key: Directed link, Value: Gurobi continuous variable
		#counted_flow = {}
		counted_F = {}	#Key: Directed link, Value: Gurobi continuous variable
		flows, counted_F = plasmids_preprocessing.flow_vars(m, links, flows, counted_F)
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
		#A link 'e' is in the plasmid only if both its endpoints are in the plasmid.
		for e in links:	
			end1, end2 = e[0], e[1]
			c1, c2 = end1[0], end2[0]
			if c1 != 'S' and c1 != 'T':
				m.addConstr(links[e] <= contigs[c1], "link_ubd-"+str(e)+'-by-'+str(c1))
				constraint_count += 1
			if c2 != 'S' and c2 != 'T':	
				m.addConstr(links[e] <= contigs[c2], "link_ubd-"+str(e)+'-by-'+str(c2))
				constraint_count += 1
					
		#Constraint type 2
		#An extremity is in the plasmid only if at least one link is incident on it.
		M = -100	#Big-M lower bound constant
		b = {}		#Auxilliary binary variable dictionary for each extremity
		for c in contigs:
			b[c] = m.addVar(vtype=GRB.BINARY, name='aux-var-'+str(c))
			x1, x2 = (c, 'h'), (c, 't')	
			expr = LinExpr()
			link_count = 0
			for link in extr_dict[x1]:
				link_count += 1
				expr.addTerms(1, links[link])
			for link in extr_dict[x2]:
				link_count += 1
				expr.addTerms(1, links[link])				
			m.addConstr(contigs[c] <= expr, "extr_ubd-"+str(c)+'-by-expr')
			m.addConstr(contigs[c] <= 1, "extr_ubd-"+str(c)+'-by-1')
			m.addConstr(contigs[c] >= expr + M*b[c], "extr_lbd")
			m.addConstr(contigs[c] >= 1 + M*(1-b[c]), "extr_lbd")
			constraint_count += 4

		#Constraint type 3
		#A contig is considered to be a ”counted” seed if it is eligible to be a seed contig
		#and is considered to be part of the solution
		seed_expr = LinExpr()
		for c in contigs:
			seed_expr.addTerms(contigs_dict[c]['Seed'], contigs[c])
		m.addConstr(seed_expr >= 1, "at-least-one-seed")
		constraint_count += 1

		#Constraint type 4
		#'F' should equal the flow out of 'S' and into 'T'. 
		#Exactly one edge exits 'S' and exactly one enters 'T'.
		xS_expr = LinExpr()
		xT_expr = LinExpr()
		fS_expr = LinExpr()
		fT_expr = LinExpr()
		#print("TESTING", len(links))
		for e in links:
			#print(e)
			if e[0] == 'S':
				xS_expr.addTerms(1, links[e])
				fS_expr.addTerms(1, flows[e])
			if e[1] == 'T':
				xT_expr.addTerms(1, links[e])
				fT_expr.addTerms(1, flows[e])
		m.addConstr(xS_expr == 1, "one-edge-from-S")
		m.addConstr(xT_expr == 1, "one-edge-into-T")
		m.addConstr(F == fS_expr, "flow-from-S")		
		m.addConstr(F == fT_expr, "flow-into-T")
		constraint_count += 4

		m.addConstr(F >= LBD_rd, "minimum flow (to avoid zero flow)")

		#Constraint types 5 and 6
		#6. Conservation constraints
		#	Flow into ('u','h') (resp. ('u','t') ) should be equal to flow out of ('u','t') (resp. ('u','h') ).
		#7. Capacity constraints
		#	The maximum flow into a vertex should be at most the capacity (read depth) of the vertex itself.
		#	The maximum flow through an edge has to be at most the capacity (capacities[e]) of the edge. 
		for c in contigs:
			ext1, ext2 = (c, 'h'), (c, 't')
			hin_flow = LinExpr()	#flow into head extremity
			hout_flow = LinExpr()	#flow out of head extremity
			tin_flow = LinExpr()	#flow into tail extremity
			tout_flow = LinExpr()	#flow out of tail extremity	

			#hin_count = LinExpr()	#flow into head extremity
			#hout_count = LinExpr()	#flow out of head extremity
			#tin_count = LinExpr()	#flow into tail extremity
			#tout_count = LinExpr()	#flow out of tail extremity	
			print("\n\n\n")
			#TODO: incoming and outgoing should be updated with graph update to avoid 'if' conditions
			if ext1 in incoming:
				for e in incoming[ext1]:
					hin_flow.addTerms(1, flows[e])
					#hin_count.addTerms(1, links[e])
					print(c, 'h', e)
				print("\n")		
			if ext1 in outgoing:
				for e in outgoing[ext1]:
					hout_flow.addTerms(1, flows[e])
					#hout_count.addTerms(1, links[e])
					print(c, 'h', e)
				print("\n")	
			if ext2 in incoming:
				for e in incoming[ext2]:
					tin_flow.addTerms(1, flows[e])
					#tin_count.addTerms(1, links[e])
					print(c, 't', e)
				print("\n")	
			if ext2 in outgoing:
				for e in outgoing[ext2]:
					tout_flow.addTerms(1, flows[e])
					#tout_count.addTerms(1, links[e])
					print(c, 't', e)
			m.addConstr(hin_flow == tout_flow, "contig-"+c+"-flow-conservation-h2t")
			m.addConstr(tin_flow == hout_flow, "contig-"+c+"-flow-conservation-t2h")

			#m.addConstr(hin_count == tout_count, "contig-"+c+"-count-conservation-h2t")
			#m.addConstr(tin_count == hout_count, "contig-"+c+"-count-conservation-t2h")			

			m.addConstr(hin_flow + tin_flow <= contigs_dict[c]['Read_depth'], "cap-"+c)	
			constraint_count += 3

			#Temporary constraint to avoid self loops
			eh2t = (ext1, ext2)
			et2h = (ext2, ext1)
			if eh2t in links:
				m.addConstr(links[eh2t] == 0, "self-loop-"+str(eh2t))
			if et2h in links:
				m.addConstr(links[et2h] == 0, "self-loop-"+str(et2h))				

		for e in links:
			m.addConstr(flows[e] <= capacities[e]*links[e], "cap-"+str(e))		
			constraint_count += 1

		#Constraint types 7 and 8
		#8. The overall flow 'F' through link 'e' is ”counted” only if 'e' is part of the solution.
		#9. The overall flow 'F' cannot exceed the flow through any active link 'e'.
		for e in links:
			m.addConstr(counted_F[e] <= UBD_rd*links[e], "xF1")
			m.addConstr(counted_F[e] <= F, "xF2")
			m.addConstr(counted_F[e] >= F - (1-links[e])*UBD_rd, "xF3")
			m.addConstr(counted_F[e] >= 0, "xF4")

			m.addConstr(counted_F[e] == flows[e], "xF5")			
			constraint_count += 5

		#Constraint type 9
		#Handling the GC-content term in the objective function
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
		constraint_count += 1

		'''
		#Constraint type 11
		#Preliminary path constraints. Note that we also account for the vertices 'S' and 'T'.
		c_expr = LinExpr()
		for c in contigs:
			c_expr.addTerms(1, contigs[c])
		l_expr = LinExpr()
		for e in links:
			l_expr.addTerms(1, links[e])
		m.addConstr(c_expr == l_expr - 1, "path_simple")					
		constraint_count += 1
		'''

		#Running the MILP
		start = time.time()
		m.optimize()
		stop = time.time()
		duration = stop - start

		#Message if solution not obtained
		if m.status == GRB.Status.INFEASIBLE:
			print ('The model cannot be solved because it is infeasible')
		elif m.status == GRB.Status.UNBOUNDED:
			print ('The model cannot be solved because it is unbounded')
		elif m.status == GRB.Status.INF_OR_UNBD:
			print ('The model cannot be solved because it is infeasible or unbounded ')

		#Storing Irreducible Inconsistent Subsystem in case solution is not obtained
		if m.status == GRB.Status.INF_OR_UNBD or m.status == GRB.Status.INFEASIBLE:
			print ('The model cannot be solved because it is infeasible')
			m.computeIIS()
			m.write("m.ilp")
			for con in m.getConstrs():
				if con.IISConstr:
					print('%s' % con.constrName) 
			exit (1)

		print("Solution:\n")
		m.printAttr('x')

		#Flow zero condition
		if F.x == 0:
			exit (1)			

		#m.write(os.path.join(output_folder, "m.sol"))

		#Retrieving plasmid bin
		output_fasta_file = open(os.path.join(output_folder, output_fasta), "a")
		ques_fasta_file = open(os.path.join(output_folder, ques_fasta), "a")
		score_file = open(os.path.join(output_folder, score_filename), "a")
		contigs_file = open(os.path.join(output_folder, output_contigs), "a")
		ques_contigs_file = open(os.path.join(output_folder, ques_contigs), "a")
		components_file = open(os.path.join(output_folder, components), "a")

		#Recording variable values (mainly for debugging purposes)
		plasmid_length = 0
	

		#for c in contigs:
		#	var_file.write(contigs[c].VarName +"\t"+ str(contigs[c].x)+"\n")
		#for e in links:
		#	var_file.write(links[e].VarName +"\t"+ str(links[e].x)+"\n")
		#	var_file.write(flows[e].VarName +"\t"+ str(flows[e].x)+"\n")
		#	var_file.write(counted_F[e].VarName +"\t"+ str(counted_F[e].x)+"\n")

		#Post-processing: Determining if plasmid bin is putative or questionable
		#Recording the plasmid bin
		for c in contigs:
			if contigs[c].x > 0:
				plasmid_length += contigs_dict[c]['Length']
		
		solution_links = set()	#Set of ordered links in the solution
		solution_seq = ''		#Sequence of contigs in plasmid bin
		from_ext = {}			#Key: Extremity, Value: List of extremities reached from this extremity by a link
		to_ext = {}				#Key: Extremity, Value: List of extremities from which this extremity is reached

		for e in links:
			if links[e].x > 0:	
				end1, end2 = e[0], e[1]		
				solution_links.add(e)
				if end1 not in from_ext:
					from_ext[end1] = [end2]
				else:
					from_ext[end1].append(end2)
				if end2 not in to_ext:
					to_ext[end2] = [end1]
				else:
					to_ext[end2].append(end1)		

		#Sorting into putative and questionable
		if plasmid_length >= 1500:
			var_file = open(os.path.join(output_folder, var_vals), "a")
			var_file.write("\nIteration "+str(n_iter)+":\n")
			for v in m.getVars():
				if v.X > 0:
					var_file.write(v.varName+"\t"+str(v.X)+"\n")

			output_fasta_file.write(">plasmid_"+str(n_iter)+"\t"+"length="+str(plasmid_length)+"\n")
			contigs_file.write("#Contig\tPrPl\tFlow\tRead_depth\tLength\n")		
			for c in contigs:
				if contigs[c].x > 0.5:			
					contigs_file.write("# "+c+"\t"+str(contigs_dict[c]['PrPl'])+"\t"+str(F.x)+"\t"+str(contigs_dict[c]['Read_depth'])+"\t"+str(contigs_dict[c]['Length'])+"\n")				
			contigs_file.write("\n")

			contigs_file.write("#Link\tCapacity\n")
			for e in links:
				if links[e].x > 0:
					contigs_file.write("# "+str(e)+"\t"+str(capacities[e])+"\n")
			contigs_file.write("\n")

			#Recording objective function scores
			GC_sum = 0
			class_sum = 0
			for c in contigs:
				gc_c_sum = 0
				for b in contig_GC[c]:
					GC_sum += -alpha2*(1-gc_probs[c][b])*contig_GC[c][b].x
					gc_c_sum += -alpha2*(1-gc_probs[c][b])*contig_GC[c][b].x
				class_sum += alpha3*contigs_dict[c]['log_ratio']*contigs[c].x
				lr_c = contigs_dict[c]['log_ratio']*contigs[c].x
				score_file.write("putative_plasmid_"+str(n_iter)+"\t\t"+str(c)+"\t"+str(gc_c_sum)+"\t"+str(lr_c)+"\n")
			score_file.write("Overall\nputative_plasmid_"+str(n_iter)+"\t\t"+str(m.objVal)+"\t"+str(F.x)+"\t"+str(GC_sum)+"\t"+str(class_sum)+"\n")
			n_iter += 1	

		#Updating assembly graph and formulation
		for c in contigs:
			if contigs[c].x > 0:
				contigs_dict[c]['Read_depth'] = max(0, contigs_dict[c]['Read_depth'] - F.x)
				if c in seeds_set and contigs_dict[c]['Read_depth'] <= 0.5:
					contigs_dict[c]['Seed'] = 0
					seeds_set.remove(c)

				if contigs_dict[c]['Read_depth'] <= 0.05:
					del(contigs_dict[c])	
					del(gc_probs[c])

					temp_list = []
					for e in links_list:						
						if (c, 'h') not in e and (c, 't') not in e:
							temp_list.append(e)
					links_list = temp_list			

		#For debugging purposes
		#seeds_set = set()				