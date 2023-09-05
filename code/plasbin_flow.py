#!/usr/bin/env python

"""
PlasBin-flow main script
"""

from re import L
from gurobipy import *
from sys import argv
import os
import time
from random import randint
import argparse
import logging
import networkx as nx

import model_setup

from data_utils import (
    LEN_KEY,
    SEED_KEY,
    COV_KEY,
    SCORE_KEY,
    UNICYCLER_TAG,
    DEFAULT_SOURCE,
    DEFAULT_SINK,
    DEFAULT_HEAD_STR,
    DEFAULT_TAIL_STR,
    DEFAULT_SEED_LEN_THRESHOLD,
    DEFAULT_SEED_SCORE_THRESHOLD,
    read_ctgs_data,
    get_seeds,
    read_gc_data,
    read_links_data,
    get_capacities,
    log_data
)
from log_errors_utils import (
    CustomException,
    process_exception
)

SOURCE = DEFAULT_SOURCE
SINK = DEFAULT_SINK

DEFAULT_SCORE_OFFSET = 0.5
DEFAULT_ALPHA1 = 1
DEFAULT_ALPHA2 = 1
DEFAULT_ALPHA3 = 1
DEFAULT_RMITER = 50
DEFAULT_MIN_PLS_LEN = 1500
DEFAULT_GUROBI_MIP_GAP = 0.05
DEFAULT_GUROBI_TIME_LIMIT = 2400

if __name__ == "__main__":
    #Parsing arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-ag", help="Path to assembly graph file")
    parser.add_argument("-gc", help="Path to GC probabilities file")
    parser.add_argument("-score", help="Path to plasmid score file")
    parser.add_argument("-out_dir", help="Path to output dir")
    parser.add_argument("-out_file", help="Name of output file")
    parser.add_argument("-log_file", help="Path to log file")
    parser.add_argument("-p", type=float, default=DEFAULT_SCORE_OFFSET, help="Offset plasmid score term")
    parser.add_argument("-alpha1", nargs='?', const=1, type=int, default=DEFAULT_ALPHA1, help="Weight of flow term")
    parser.add_argument("-alpha2", nargs='?', const=1, type=int, default=DEFAULT_ALPHA2, help="Weight of GC content term")
    parser.add_argument("-alpha3", nargs='?', const=1, type=int, default=DEFAULT_ALPHA3, help="Weight of log probabilities term")	
    parser.add_argument("-rmiter", nargs='?', const=1, type=int, default=DEFAULT_RMITER, help="Number of iterations to remove circular components")
    parser.add_argument("-assembler", nargs='?', const=1, type=str, default=UNICYCLER_TAG, help="Name of assembler (unicycler/skesa)")
    parser.add_argument("-seed_len", nargs='?', const=1, type=int, default=DEFAULT_SEED_LEN_THRESHOLD, help="Seed length threshold")
    parser.add_argument("-seed_score", nargs='?', const=1, type=float, default=DEFAULT_SEED_SCORE_THRESHOLD, help="Seed plasmid score threshold")
    parser.add_argument("-gc_intervals", nargs='?', const=1, default=None, help="GC intervals file")
    parser.add_argument("-min_pls_len", nargs='?', const=1, type=int, default=DEFAULT_MIN_PLS_LEN, help="Minimum plasmid length to be reported")
    parser.add_argument("-gurobi_mip_gap", nargs='?', const=1, type=float, default=DEFAULT_GUROBI_MIP_GAP, help="MIPGap parameter for Gurobi")
    parser.add_argument("-gurobi_time_limit", nargs='?', const=1, type=int, default=DEFAULT_GUROBI_TIME_LIMIT, help="Time limit for Gurobi (in seconds)")        
    
    args = parser.parse_args()
    
    output_dir = args.out_dir
    output_file = args.out_file
    log_file = args.log_file
    assembly_file = args.ag
    gc_prob_file = args.gc
    score_file = args.score
    p = float(args.p)
    alpha1 = float(args.alpha1)
    alpha2 = float(args.alpha2)
    alpha3 = float(args.alpha3)
    rmiter = args.rmiter
    assembler = args.assembler
    seed_len = args.seed_len
    seed_score = args.seed_score
    gc_int_file = args.gc_intervals
    min_pls_len = args.min_pls_len
    gurobi_mip_gap = args.gurobi_mip_gap
    gutobi_time_limit = args.gurobi_time_limit
    
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s'
    )

    #Naming and creating output files
    #output_folder = output_dir + '/' + ratios
    output_folder = output_dir
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

    # Reading data
    contigs_dict = read_ctgs_data(
        assembly_file, score_file,
        assembler=assembler, gfa_gzipped=True
    )
    seeds_set = get_seeds(contigs_dict, seed_len=seed_len, seed_score=seed_score)
    gc_probs, gc_pens = read_gc_data(gc_prob_file, gc_int_file)
    links_list = read_links_data(assembly_file, gfa_gzipped=True)
    try:
        gc_prob_ctgs_list = sorted(gc_probs.keys())
	ctgs_data_ctgs_list = sorted(ctgs_data_dict.keys())
	if gc_prob_ctgs_list != ctgs_data_ctgs_list:
	    raise CustomException(f'Inconsistent contigs sets')
    except Exception as e:
	process_exception(f'Files {assembly_file} and {gc_prob_file}')
    log_data(contigs_dict, links_list, assembly_file, score_file)
    
    # contigs_dict = {}   #Key: contig IDs, Values: Contig attributes (provided as or derived from input)
    # links_list = []     #List of unordered links: Each link is a pair of extrmeities (e.g. (('1',DEFAULT_HEAD_STR),('2',DEFAULT_TAIL_STR)))
    # gc_probs = {}       #Key: contigs IDs, Values: Probability for each GC bin
    # gc_pens = {}		#Key: contigs IDs, Values: Penalty for each GC bin
    # capacities = {}     #Key: 
    # contigs_dict, links_list = get_data.get_ag_details(assembly_file, contigs_dict, links_list)
    # contigs_dict = get_data.get_gene_coverage(mapping_file, contigs_dict)
    # gc_probs, gc_pens = get_data.get_gc_probs(gc_file, gc_probs, gc_pens)

    n_comp = 0
    while len(seeds_set) > 0:
        #For consistency, both extremities for a contig should be part of exactly the same number of links in a plasmid.
        #For each extremity, we make a list of links that involve the extremity.
        extr_dict = {}	#Key: Extremity (e.g. ('1',DEFAULT_HEAD_STR)), Value: List of unordered links incident on extremity
	#List of incoming and outgoing edges for each extremity
	#Including links from S and links to T
        incoming, outgoing = {}, {}   #Key: Extremity, Value: List of ordered links in/out of extremity

        UBD_rd, LBD_rd = 0, 100
        for c in contigs_dict:
            UBD_rd = max(UBD_rd, contigs_dict[c][COV_KEY])
            LBD_rd = min(LBD_rd, contigs_dict[c][COV_KEY])

            ext1, ext2 = (c, DEFAULT_HEAD_STR), (c, DEFAULT_TAIL_STR)
            extr_dict[ext1], extr_dict[ext2] = [], []
            incoming[ext1], incoming[ext2] = [], []
            outgoing[ext1], outgoing[ext2] = [], []
            if c in seeds_set:
                extr_dict[ext1], extr_dict[ext2] = [(SOURCE,ext1)], [(SOURCE,ext2)]
                incoming[ext1], incoming[ext2] = [(SOURCE, ext1)], [(SOURCE, ext2)]
                extr_dict[ext1], extr_dict[ext2] = [(ext1,SINK)], [(ext2,SINK)]
            outgoing[ext1], outgoing[ext2] = [(ext1, SINK)], [(ext2, SINK)]

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

            capacities = get_capacities(links_list, contigs_dict)

            #-----------------------------------------------
	    #Initializing the ILP
            m = Model("Plasmids")
            m.params.LogFile= os.path.join(output_folder,'m.log')
            m.setParam(GRB.Param.TimeLimit, gurobi_time_limit)
            m.setParam(GRB.Param.MIPGap, gurobi_mip_gap)

	    #Initializing variables
            contigs = {}	#Key: Contig (e.g. '1'), Value: Gurobi binary variable
            contigs = model_setup.contig_vars(m, contigs_dict, contigs)

            links = {}		#Key: Directed link from one extremity to another
            #(e.g. (('1',DEFAULT_HEAD_STR),('2',DEFAULT_TAIL_STR)) ), Value: Gurobi binary variable
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
            # alpha1, alpha2, alpha3 = float(alpha1), float(alpha2), float(alpha3)

            expr.addTerms(alpha1, F)
            for c in contigs:
                for b in plas_GC:
                    expr.addTerms(alpha2*(gc_pens[c][b]), contig_GC[c][b])
                expr.addTerms(alpha3*(contigs_dict[c][SCORE_KEY] - p), contigs[c])
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
	    #'F' should equal the flow out of SOURCE and into SINK. 
	    #Exactly one edge exits SOURCE and exactly one enters SINK.
            m, constraint_count = model_setup.min_flow_constr(m, links, flows, F, LBD_rd, constraint_count)

	    #Constraint types 5 and 6
	    #6. Conservation constraints
	    #	Flow into ('u',DEFAULT_HEAD_STR) (resp. ('u',DEFAULT_TAIL_STR) ) should be equal to flow out of ('u',DEFAULT_TAIL_STR) (resp. ('u',DEFAULT_HEAD_STR) ).
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
            while extra_comps >= 1 and iter_count <= rmiter:
		#Running the MILP
                start = time.time()
                m.optimize()
                stop = time.time()
                duration = stop - start
                logging.info(f'MILP\tIteration {iter_count}: {duration}')

                iter_count += 1

		#Message if solution not obtained
                if m.status == GRB.Status.INFEASIBLE:
                    logging.warning(f'MILP\tThe model cannot be solved because it is infeasible')
                elif m.status == GRB.Status.UNBOUNDED:
                    logging.warning(f'MILP\tThe model cannot be solved because it is unbounded')
                elif m.status == GRB.Status.INF_OR_UNBD:
                    logging.warning(f'MILP\tThe model cannot be solved because it is infeasible or unbounded ')

		#Storing Irreducible Inconsistent Subsystem in case solution is not obtained
                if m.status == GRB.Status.INF_OR_UNBD or m.status == GRB.Status.INFEASIBLE:
                    print(f'MILP\tStoring Irreducible Inconsistent Subsystem in m.ilp')
                    m.computeIIS()
                    m.write("m.ilp")
                    for con in m.getConstrs():
                        if con.IISConstr:
                            print(f'{con.constrName}')
                    exit(1)

                print("Solution:\n")
                m.printAttr('x')

		#Flow zero condition
                if F.x == 0:
                    exit(1)			

		#Finding components in the solution
                G = nx.DiGraph()
                for c in contigs:
                    if contigs[c].x > 0:
                        G.add_node(c)
                    G.add_node(SOURCE)
                    G.add_node(SINK)

                for e in links:
                    if links[e].x > 0:
                        end1, end2 = e[0], e[1]
                        if end1 == SOURCE:
                            c1 = SOURCE
                            ext1 = None
                        else:
                            c1 = end1[0]
                            ext1 = end1[1]
                        if end2 == SINK:
                            c2 = SINK
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
                        if node != SOURCE and node != SINK:
                            comp_len += contigs_dict[node][LEN_KEY]
                    if SOURCE in comp:
                        ST_comp = G.subgraph(comp)
                        print("Edges:",ST_comp.edges())
                        for edge in ST_comp.edges:
                            print(nx.get_edge_attributes(G,'extremities')[edge])
                    else:
                        disconn_comp = G.subgraph(comp)
                        #Muting individual edges
                        for edge in disconn_comp.edges:
                            exts = nx.get_edge_attributes(G,'extremities')[edge]
                            e = ((edge[0],exts[0]),(edge[1],exts[1]))
                            m.addConstr(links[e] == 0, "muted_edge-"+str(e))	

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

	    #Recording the plasmid bin if the plasmid is long enough
            for c in contigs:
                if contigs[c].x > 0:
                    plasmid_length += contigs_dict[c][LEN_KEY]

            if plasmid_length >= min_pls_len:
		#Recording objective function scores
                GC_sum = 0
                gd_sum = 0
                for c in contigs:
                    gc_c_sum = 0
                    for b in contig_GC[c]:
                        GC_sum += alpha2*(gc_pens[c][b])*contig_GC[c][b].x
                        gc_c_sum += alpha2*(gc_pens[c][b])*contig_GC[c][b].x
                    gd_sum += alpha3*(contigs_dict[c][SCORE_KEY]-0.5)*contigs[c].x
                    lr_gd = (contigs_dict[c][SCORE_KEY]-0.5)*contigs[c].x

		#Recording components in the solution
                G = nx.DiGraph()
                for c in contigs:
                    if contigs[c].x > 0:
                        G.add_node(c)
                G.add_node(SOURCE)
                G.add_node(SINK)

                for e in links:
                    if links[e].x > 0:
                        end1, end2 = e[0], e[1]
                        if end1 == SOURCE:
                            c1 = SOURCE
                            ext1 = None
                        else:
                            c1 = end1[0]
                            ext1 = end1[1]
                        if end2 == SINK:
                            c2 = SINK
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
                    comp_len = 0
                    for node in comp:
                        if node != SOURCE and node != SINK:
                            if node not in pbf_bins[n_iter]['Contigs']:
                                pbf_bins[n_iter]['Contigs'][node] = 0

	    #Updating assembly graph and formulation
            for e in flows:
                if e[1] != SINK:
                    c = e[1][0]
                    contigs_dict[c][COV_KEY] = max(0, contigs_dict[c][COV_KEY] - flows[e].x)
                    if c in pbf_bins[n_iter]['Contigs']:
                        pbf_bins[n_iter]['Contigs'][c] += round(flows[e].x / F.x, 2)

            for c in contigs:
                if contigs[c].x > 0:
                    if c in seeds_set and contigs_dict[c][COV_KEY] <= 0.5:
                        contigs_dict[c][SEED_KEY] = 0
                        seeds_set.remove(c)

                    if contigs_dict[c][COV_KEY] <= 0.05:
                        del(contigs_dict[c])
                        del(gc_probs[c])
                        del(gc_pens[c])

                        temp_list = []
                        for e in links_list:
                            if (c, DEFAULT_HEAD_STR) not in e and (c, DEFAULT_TAIL_STR) not in e:
                                temp_list.append(e)
                        links_list = temp_list	


    output_bins.write('#Pls_ID\tFlow\tGC_bin\t\tContigs\n')
    for p in pbf_bins:
        fval = "%.2f" %pbf_bins[p]['Flow']
        gcb = pbf_bins[p]['GC_bin']
        print(p, fval, gcb, pbf_bins[p]['Contigs'])

        output_bins.write('P'+str(p)+'\t\t'+str(fval)+'\t'+str(gcb)+'\t\t')
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
