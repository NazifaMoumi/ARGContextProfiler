#!/usr/bin/env python

import os
import json
import networkx as nx
from src import config
from src import graph_utils, diamond_filter, annotation, path_extraction, context_extraction
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run ARG Context Pipeline")
    
    # File paths and directories
    parser.add_argument("--graph", default=config.GRAPH_FILE_PATH, help="Path to the assembly graph in .fastg format")
    parser.add_argument("--db", default=config.DB_PATH, help="Path to the reference DB diamond file")
    parser.add_argument("--db_fasta", default=config.DB_PATH_FASTA, help="Path to the reference DB FASTA file")
    parser.add_argument("--output", default=config.OUTPUT_PATH, help="Directory for output results")
    
    # Parameters
    parser.add_argument("--overlap", type=int, default=config.OVERLAP, help="Overlap length between nodes")
    parser.add_argument("--radius", type=int, default=config.RADIUS, help="Radius for context extraction")
    parser.add_argument("--gene_cov", type=float, default=config.GENE_COV, help="Minimum gene coverage ratio")
    parser.add_argument("--max_gap_length", type=int, default=config.MAX_GAP_LENGTH, help="Maximum allowed gap length")
    parser.add_argument("--max_node_dist", type=int, default=config.MAX_NODE_DIST, help="Maximum node distance")
    parser.add_argument("--max_gap_node", type=int, default=config.MAX_GAP_NODE, help="Maximum allowed gap nodes")
    parser.add_argument("--min_context_len", type=int, default=config.MIN_CONTEXT_LEN, help="Minimum context length")
    parser.add_argument("--group_cov_ratio", type=float, default=config.GROUP_COV_RATIO, help="Group coverage ratio")
    parser.add_argument("--assembler", default=config.ASSEMBLER, help="Assembler used (e.g., 'metaspades' or 'megahit')")
    
    return parser.parse_args()

def main():

    args = parse_arguments()
    
    graph_file = args.graph
    db_path = args.db
    db_path_fasta = args.db_fasta
    output_path = args.output
    overlap = args.overlap
    radius = args.radius
    gene_cov = args.gene_cov
    max_gap_length = args.max_gap_length
    max_node_dist = args.max_node_dist
    max_gap_node = args.max_gap_node
    min_context_len = args.min_context_len
    group_cov_ratio = args.group_cov_ratio
    assembler = args.assembler

    print('Input graph:', graph_file)
    print('Output path:', output_path)
    
    # (Optional) Format graph if needed.
    # if assembler == 'megahit':
    #     formatted_graph = graph_utils.format_graph(graph_file)
    # elif assembler == 'metaspades':
    #     formatted_graph = graph_file  # or use formatted version
    # else:
    #     print('assembler not recognized, terminating...')
    #     return

    # # Create sequence file from graph.
    # seq_dict = graph_utils.create_seq_file_from_graph(formatted_graph, output_path)
    
    # # Create graph from sequences.
    # g = graph_utils.create_graph(seq_dict, formatted_graph)
    # print(f"Graph built with {g.number_of_nodes()} nodes and {g.number_of_edges()} edges")
    
    # # Save graph if needed.
    # nx.write_gexf(g, os.path.join(output_path, 'nx_graph.gexf'))
    g = nx.read_gexf(output_path + 'nx_graph.gexf')
    print('graph read with ' + str(g.number_of_nodes()) + ' nodes and ' + str(g.number_of_edges()) + ' edges')

    # create_filtered_seq_file(g1)

    # Annotate ARGs.
    ARG_dict, ARG_target_dict = annotation.annotate_ARG(db_path, output_path)
    print("Total ARGs:", len(ARG_dict))
    with open(os.path.join(output_path, 'ARG_dict_content.json'), 'w') as f:
        json.dump(ARG_dict, f, indent=4)
    
    ARG_dict = dict(sorted(ARG_dict.items(), key=lambda item: len(item[1])))

    # Extract ARG paths.
    ARG_path_dict, ARG_path_list = path_extraction.extract_ARG_paths_seq_gap(
        g, ARG_dict, db_path_fasta, overlap, gene_cov, max_gap_length, max_gap_node, group_cov_ratio, max_node_dist)
    
    # Optionally remove previous result files.
    for fname in ['left_seq.fasta', 'right_seq.fasta', 'contexts.txt', 'ARG_seq.fasta']:
        fpath = os.path.join(output_path, f"result_contexts_cov_{gene_cov}_radius_{radius}", fname)
        if os.path.exists(fpath):
            os.remove(fpath)
    
    # Extract genomic contexts.
    context_extraction.extract_ARG_contexts(g, ARG_dict, ARG_path_dict, radius, overlap, output_path, gene_cov, min_context_len)
    print("Gene context extraction done.")

if __name__ == '__main__':
    main()
