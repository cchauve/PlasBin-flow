__author__ = 'amane'

#-------------------
#USAGE: 
#time python plasbin_flow.py --ag assembly.gfa 
# 				--gc gc_probs.csv --map gene_contig_mapping.csv \
#				--outdir output_dir --outfile output_file \
# 				--alpha1 alpha_1 --alpha2 alpha_2 --alpha3 alpha_3 \
# 				--offset offset --rmiter rmiter --unique unique

from re import L
from gurobipy import *
from sys import argv
import os
import time
from random import randint
import argparse
#import preprocessing
import get_data
import model_setup
import networkx as nx

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
	parser.add_argument("--map", help="Path to gene to contig mapping file")
	parser.add_argument("--p", type=float, default = 0.5, help="Offset gd / plasmid score term")
	parser.add_argument("--outdir", help="Path to output dir")
	parser.add_argument("--outfile", help="Name of output file")
	parser.add_argument("--alpha1", nargs='?', const = 1, type=int, default = 1, help="Weight of flow term")
	parser.add_argument("--alpha2", nargs='?', const = 1, type=int, default = 1, help="Weight of GC content term")
	parser.add_argument("--alpha3", nargs='?', const = 1, type=int, default = 1, help="Weight of log probabilities term")	
	parser.add_argument("--rmiter", nargs='?', const = 1, type=int, default = 50, help="Number of iterations to remove circular components")
	parser.add_argument("--unique", nargs='?', const = 1, type=int, default = 1, help="Unique bins")

	args = parser.parse_args()

	output_dir = args.outdir
	output_file = args.outfile
	assembly_file = args.ag
	gc_file = args.gc
	mapping_file = args.map
	p = float(args.p)
	alpha1 = args.alpha1
	alpha2 = args.alpha2
	alpha3 = args.alpha3
	rmiter = int(args.rmiter)
	unique = int(args.unique)

	#Naming and creating output files
	ratios = str(alpha1) + '.' + str(alpha2) + '.' + str(alpha3)
	output_folder = output_dir + '/' + ratios
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	output_bins = open(os.path.join(output_folder, output_file), "w")

	
	
	#output_fasta = 'putative_plasmids.fasta'
	#score_filename = 'MILP_objective.csv'
	#var_vals = 'var_values.txt'

	#Dictionary for storing bins
	#Key = Plasmid id, Value = Dictionary of attributes
	#Flow (float), GC bin id, List of contigs with multiplicities
	pbf_bins = {}
	n_iter = 0

	#-----------------------------------------------
	#Main program
	contigs_dict = {}   #Key: contig IDs, Values: Contig attributes (provided as or derived from input)
	links_list = []     #List of unordered links: Each link is a pair of extrmeities (e.g. (('1','h'),('2','t')))
	seeds_set = set()   #Set of contig IDs from seeds file
	gc_probs = {}       #Key: contigs IDs, Values: Probability for each GC bin
	gc_pens = {}		#Key: contigs IDs, Values: Penalty for each GC bin
	capacities = {}     #Key: 

	contigs_dict, links_list = get_data.get_ag_details(assembly_file, contigs_dict, links_list)
	contigs_dict = get_data.get_gene_coverage(mapping_file, contigs_dict)
	
	seeds_set = get_data.get_seeds(contigs_dict, seeds_set)
	gc_probs, gc_pens = get_data.get_gc_probs(gc_file, gc_probs, gc_pens)

	for c in contigs_dict:  #Assigning seeds
		if c in seeds_set:
			contigs_dict[c]['Seed'] = 1
		else:
			contigs_dict[c]['Seed'] = 0

	#Keeping a log of details about contigs in the assembly graph
	'''
	input_details = 'details.csv'
	details_file = open(os.path.join(output_folder, input_details), "w")
	details_file.write("Contig"+"\t"+"Read_depth"+ "\t"+"GC_cont"+"\t"+ "Length"+"\t"+"Density"+"\n")
	for c in contigs_dict:
		details_file.write(c+"\t"+str(contigs_dict[c]['Read_depth'])+ "\t"+str(contigs_dict[c]['GC_cont'])+"\t"+ str(contigs_dict[c]['Length'])+"\t"+str(contigs_dict[c]['Gene_coverage'])+"\n")
		details_file.write("Link list:"+"\t")
		for link in links_list:
			if c == link[0][0] or c == link[1][0]:
				details_file.write(str(link)+"\t")
		details_file.write("\n")
	'''

	#output_fasta_file = open(os.path.join(output_folder, output_fasta), "w")
	#score_file = open(os.path.join(output_folder, score_filename), "w")
	#contigs_file = open(os.path.join(output_folder, output_contigs), "w")
	#var_file = open(os.path.join(output_folder, var_vals), "w")

	
	n_comp = 0

	while len(seeds_set) > 0:
		#For consistency, both extremities for a contig should be part of exactly the same number of links in a plasmid.
		#For each extremity, we make a list of links that involve the extremity.
		extr_dict = {}	#Key: Extremity (e.g. ('1','h')), Value: List of unordered links incident on extremity
		#List of incoming and outgoing edges for each extremity
		#Including links from S and links to T
		incoming, outgoing = {}, {}	#Key: Extremity, Value: List of ordered links in/out of extremity

		UBD_rd, LBD_rd = 0, 100
		for c in contigs_dict:
			UBD_rd = max(UBD_rd, contigs_dict[c]['Read_depth'])
			LBD_rd = min(LBD_rd, contigs_dict[c]['Read_depth'])

			ext1, ext2 = (c, 'h'), (c, 't')
			extr_dict[ext1], extr_dict[ext2] = [], []
			incoming[ext1], incoming[ext2] = [], []
			outgoing[ext1], outgoing[ext2] = [], []
			if c in seeds_set:
				extr_dict[ext1], extr_dict[ext2] = [('S',ext1)], [('S',ext2)]
				incoming[ext1], incoming[ext2] = [('S', ext1)], [('S', ext2)]
			extr_dict[ext1], extr_dict[ext2] = [(ext1,'T')], [(ext2,'T')]
			outgoing[ext1], outgoing[ext2] = [(ext1, 'T')], [(ext2, 'T')]

			for link in links_list:
				if link[0] == ext1 or link[1] == ext1:
					extr_dict[ext1].append(link) 
					extr_dict[ext1].append(link[::-1]) 
				if link[0] == ext2 or link[1] == ext2:
					extr_dict[ext2].append(link)
					extr_dict[ext2].append(link[::-1])

		for link in links_list:
			ext1, ext2 = link[0], link[1]
			incoming[ext1].append(link[::-1])
			outgoing[ext1].append(link)
			incoming[ext2].append(link)
			outgoing[ext2].append(link[::-1])		
		
		capacities = get_data.get_caps(links_list, contigs_dict, capacities)

		print("\n\n\n\n\n")
		#-----------------------------------------------
		#Initializing the ILP
		m = Model("Plasmids")
		m.params.LogFile= os.path.join(output_folder,'m.log')
		m.setParam(GRB.Param.TimeLimit, 2400)
		m.setParam(GRB.Param.MIPGap, 0.05)

		#Initializing variables
		contigs = {}	#Key: Contig (e.g. '1'), Value: Gurobi binary variable
		contigs = model_setup.contig_vars(m, contigs_dict, contigs)

		links = {}		#Key: Directed link from one extremity to another 
						#(e.g. (('1','h'),('2','t')) ), Value: Gurobi binary variable
		links = model_setup.link_vars(m, links_list, links, contigs)

		plas_GC = {}	#Key: GC bin, Value: Gurobi binary variable
		contig_GC = {}	#Nested Dict: Key: Contig, Value: Dict of GC bins (Key: GC bin, Value: Gurobi binary variable)
		plas_GC, contig_GC = model_setup.GC_vars(m, gc_probs, plas_GC, contig_GC)

		flows = {}		#Key: Directed link, Value: Gurobi continuous variable
		counted_F = {}	#Key: Directed link, Value: Gurobi continuous variable
		flows, counted_F = model_setup.flow_vars(m, links, flows, counted_F)
		F = m.addVar(vtype=GRB.CONTINUOUS, name='overall-flow')	

		#-----------------------------------------------
		#Setting up the expression for the objective function
		expr = LinExpr()
		alpha1, alpha2, alpha3 = float(alpha1), float(alpha2), float(alpha3)

		expr.addTerms(alpha1, F)
		for c in contigs:
			for b in plas_GC:
				#expr.addTerms(-alpha2*(1-gc_probs[c][b]), contig_GC[c][b])
				expr.addTerms(alpha2*(gc_pens[c][b]), contig_GC[c][b])
			#expr.addTerms(alpha3*contigs_dict[c]['log_ratio'], contigs[c])
			expr.addTerms(alpha3*(contigs_dict[c]['Gene_coverage'] - p), contigs[c])
		m.setObjective(expr, GRB.MAXIMIZE)

		#-----------------------------------------------
		#Setting up constraints

		constraint_count = 0
		#Constraint type 1
		#A link 'e' is in the plasmid only if both its endpoints are in the plasmid.
		m, constraint_count = model_setup.link_inclusion_constr(m, links, contigs, constraint_count)
					
		#Constraint type 2
		#An extremity is in the plasmid only if at least one link is incident on it.
		m, constraint_count = model_setup.extr_inclusion_constr(m, links, contigs, extr_dict, constraint_count)

		#Constraint type 3
		#A contig is considered to be a ”counted” seed if it is eligible to be a seed contig
		#and is considered to be part of the solution
		m, constraint_count = model_setup.seed_inclusion_constr(m, contigs, contigs_dict, constraint_count)

		#Constraint type 4
		#'F' should equal the flow out of 'S' and into 'T'. 
		#Exactly one edge exits 'S' and exactly one enters 'T'.
		m, constraint_count = model_setup.min_flow_constr(m, links, flows, F, LBD_rd, constraint_count)

		#Constraint types 5 and 6
		#6. Conservation constraints
		#	Flow into ('u','h') (resp. ('u','t') ) should be equal to flow out of ('u','t') (resp. ('u','h') ).
		#7. Capacity constraints
		#	The maximum flow into a vertex should be at most the capacity (read depth) of the vertex itself.
		#	The maximum flow through an edge has to be at most the capacity (capacities[e]) of the edge. 
		m, constraint_count = model_setup.flow_conservation_constraints(m, links, contigs, flows, incoming, outgoing, capacities, contigs_dict, constraint_count)

		#Constraint types 7 and 8
		#8. The overall flow 'F' through link 'e' is ”counted” only if 'e' is part of the solution.
		#9. The overall flow 'F' cannot exceed the flow through any active link 'e'.
		m, constraint_count = model_setup.counted_flow_constr(m, links, flows, counted_F, F, UBD_rd, constraint_count)

		#Constraint type 9
		#Handling the GC-content term in the objective function
		m, constraint_count = model_setup.GC_constr(m, contig_GC, plas_GC, contigs, constraint_count)

		extra_comps = 1	#default
		iter_count = 0
		dc_count = 0
		dc_dict = {}
		while extra_comps >= 1 and iter_count <= 50:
			#Running the MILP
			start = time.time()
			m.optimize()
			stop = time.time()
			duration = stop - start
			print(duration)

			iter_count += 1

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

			#Finding components in the solution
			G = nx.DiGraph()
			for c in contigs:
				if contigs[c].x > 0:
					G.add_node(c)
			G.add_node('S')
			G.add_node('T')

			for e in links:
				if links[e].x > 0:
					end1, end2 = e[0], e[1]
					if end1 == 'S':
						c1 = 'S'
						ext1 = None
					else:
						c1 = end1[0]
						ext1 = end1[1]
					if end2 == 'T':
						c2 = 'T'
						ext2 = None
					else:
						c2 = end2[0]
						ext2 = end2[1]
					G.add_edge(c1,c2)
					nx.set_edge_attributes(G, {(c1, c2): {"extremities": (ext1, ext2)}})
			conn_comps = nx.weakly_connected_components(G)

			#components_file = open(os.path.join(output_folder, components), "a")
			comp_count = 0
			for comp in conn_comps:
				comp_count += 1
				comp_len = 0
				for node in comp:
					if node != 'S' and node != 'T':
						comp_len += contigs_dict[node]['Length']
				#print("Conn comp:", comp_count, comp_len)
				if 'S' in comp:
					ST_comp = G.subgraph(comp)
					print("Edges:",ST_comp.edges())
					for edge in ST_comp.edges:
						print(nx.get_edge_attributes(G,'extremities')[edge])
				else:
					disconn_comp = G.subgraph(comp)
					#Muting individual edges
					for edge in disconn_comp.edges:
						#print(nx.get_edge_attributes(G,'extremities')[edge])	
						exts = nx.get_edge_attributes(G,'extremities')[edge]
						e = ((edge[0],exts[0]),(edge[1],exts[1]))
						#expr.addTerms(1, links[e])
						
						m.addConstr(links[e] == 0, "muted_edge-"+str(e))	
					'''
					dc_count += 1
					dc_dict[dc_count] = m.addVar(vtype=GRB.BINARY, name='dc-binary-'+str(dc_count))
					eps = 0.001
					bigM = 1000
					expr = LinExpr()
					
					#print("Disconnected comp:",disconn_comp.edges())
					nodes = disconn_comp.nodes()
					#edges = disconn_comp.edges
					edges_into_comp = LinExpr()
					for node in nodes:
						#print(node)
						h_ext = (node, 'h')
						t_ext = (node, 't')

						for e in incoming[h_ext]:
							if e in links:
								edges_into_comp.addTerms(1, links[e])
						for e in incoming[t_ext]:		
							if e in links:
								edges_into_comp.addTerms(1, links[e])

					m.addConstr(edges_into_comp >= 0 + eps - bigM * dc_dict[dc_count], name="bigM-"+str(dc_count))			


					#m.addConstr(x >= y + eps - M * (1 - b), name="bigM_constr1")
					#m.addConstr(x <= y + M * b, name="bigM_constr2")					
					'''


			#Condition to stop iterating. 
			#If number of connected components is 1, there are no extra components, thus breaking the while loop.	
			extra_comps = comp_count - 1			

		#Plasmid bin obtained for the iteration.
		#Out of while loop to remove extra circular components. 
		#Proceeding to output

		#Retrieving plasmid bin
		#output_fasta_file = open(os.path.join(output_folder, output_fasta), "a")
		#score_file = open(os.path.join(output_folder, score_filename), "a")
		#contigs_file = open(os.path.join(output_folder, output_contigs), "a")
		#components_file = open(os.path.join(output_folder, components), "a")


		#Recording variable values (mainly for debugging purposes)
		plasmid_length = 0

		#Post-processing: Determining if plasmid bin is putative or questionable
		#Recording the plasmid bin
		for c in contigs:
			if contigs[c].x > 0:
				plasmid_length += contigs_dict[c]['Length']
		
			

		#Sorting into putative and questionable
		if plasmid_length >= 1500:

			#var_file = open(os.path.join(output_folder, var_vals), "a")
			#var_file.write("\nIteration "+str(n_iter)+":\n")
			#for v in m.getVars():
			#	if v.X > 0:
			#		var_file.write(v.varName+"\t"+str(v.X)+"\n")

			#output_fasta_file.write(">plasmid_"+str(n_iter)+"\t"+"length="+str(plasmid_length)+"\n")
			#contigs_file.write("#Contig\tDensity\tFlow\tRead_depth\tLength\n")		
			#for c in contigs:
			#	if contigs[c].x > 0.5:			
			#		contigs_file.write("# "+c+"\t"+str(contigs_dict[c]['Gene_coverage'])+"\t"+str(F.x)+"\t"+str(contigs_dict[c]['Read_depth'])+"\t"+str(contigs_dict[c]['Length'])+"\n")				
			#contigs_file.write("\n")

			#contigs_file.write("#Link\tCapacity\n")
			#for e in links:
			#	if links[e].x > 0:
			#		contigs_file.write("# "+str(e)+"\t"+str(capacities[e])+"\n")
			#contigs_file.write("\n")

			#Recording objective function scores
			GC_sum = 0
			gd_sum = 0
			for c in contigs:
				gc_c_sum = 0
				for b in contig_GC[c]:
					#GC_sum += -alpha2*(1-gc_probs[c][b])*contig_GC[c][b].x
					GC_sum += alpha2*(gc_pens[c][b])*contig_GC[c][b].x
					gc_c_sum += alpha2*(gc_pens[c][b])*contig_GC[c][b].x
				gd_sum += alpha3*(contigs_dict[c]['Gene_coverage']-0.5)*contigs[c].x
				lr_gd = (contigs_dict[c]['Gene_coverage']-0.5)*contigs[c].x
				#score_file.write("putative_plasmid_"+str(n_iter)+"\t\t"+str(c)+"\t"+str(gc_c_sum)+"\t"+str(lr_gd)+"\n")
			#score_file.write("Overall\nputative_plasmid_"+str(n_iter)+"\t\t"+str(m.objVal)+"\t"+str(F.x)+"\t"+str(GC_sum)+"\t"+str(gd_sum)+"\n\n")
				


			#Recording components in the solution
			G = nx.DiGraph()
			for c in contigs:
				if contigs[c].x > 0:
					G.add_node(c)
			G.add_node('S')
			G.add_node('T')

			for e in links:
				if links[e].x > 0:
					end1, end2 = e[0], e[1]
					if end1 == 'S':
						c1 = 'S'
						ext1 = None
					else:
						c1 = end1[0]
						ext1 = end1[1]
					if end2 == 'T':
						c2 = 'T'
						ext2 = None
					else:
						c2 = end2[0]
						ext2 = end2[1]
					#print(c1, c2)	
					G.add_edge(c1,c2)
					nx.set_edge_attributes(G, {(c1, c2): {"extremities": (ext1, ext2)}})
			conn_comps = nx.weakly_connected_components(G)

			#components_file = open(os.path.join(output_folder, components), "a")
			GC_bin = 0
			for b in plas_GC:
				if plas_GC[b].x == 1:
					GC_bin = b

			comp_count = 0
			for comp in conn_comps:
				n_iter += 1	
				pbf_bins[n_iter] = {'Flow': F.x, 'GC_bin': GC_bin, 'Contigs': {}}

				comp_count += 1
				#print(comp)
				comp_len = 0
				#components_file.write("Component "+str(comp_count)+":\t")
				for node in comp:
					if node != 'S' and node != 'T':
						if node not in pbf_bins[n_iter]['Contigs']:
							pbf_bins[n_iter]['Contigs'][node] = 0
					#components_file.write(str(node)+",")
				#components_file.write("\n")				

		#Updating assembly graph and formulation
		for e in flows:
			if e[1] != 'T':
				c = e[1][0]
				contigs_dict[c]['Read_depth'] = max(0, contigs_dict[c]['Read_depth'] - flows[e].x)
				if c in pbf_bins[n_iter]['Contigs']:
					pbf_bins[n_iter]['Contigs'][c] += round(flows[e].x / F.x, 2)

			
		for c in contigs:
			if contigs[c].x > 0:
				if c in seeds_set and contigs_dict[c]['Read_depth'] <= 0.5:
					contigs_dict[c]['Seed'] = 0
					seeds_set.remove(c)

				if contigs_dict[c]['Read_depth'] <= 0.05:
					del(contigs_dict[c])	
					del(gc_probs[c])
					del(gc_pens[c])

					temp_list = []
					for e in links_list:						
						if (c, 'h') not in e and (c, 't') not in e:
							temp_list.append(e)
					links_list = temp_list	
		

				
		#For debugging purposes
		#seeds_set = set()	
	output_bins.write('#Pls_ID\tFlow\tGC_bin\t\tContigs\n')		
	for p in pbf_bins:
		fval = "%.2f" %pbf_bins[p]['Flow']
		gcb = pbf_bins[p]['GC_bin']
		print(p, fval, gcb, pbf_bins[p]['Contigs'])
		
		output_bins.write('P'+str(p)+'\t\t'+str(fval)+'\t'+str(gcb)+'\t\t')
		#nctg = len(pbf_bins[p]['Contigs'])
		nctg = 0
		for c in pbf_bins[p]['Contigs']:
			ctg_mul = pbf_bins[p]['Contigs'][c]
			if nctg == 0:
				output_bins.write(c+':'+str(ctg_mul))
			else:
				output_bins.write(','+c+':'+str(ctg_mul))
			nctg += 1
		output_bins.write('\n')
		
	print("Out of while loop")


	#print(pbf_bins)
		