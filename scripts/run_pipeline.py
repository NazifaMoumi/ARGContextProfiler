#!/usr/bin/env python

import os
import json
import networkx as nx
from src import config
from src import graph_utils, diamond_filter, annotation, path_extraction, context_extraction

def main():
    print('Input graph:', config.GRAPH_FILE_PATH)
    print('Output path:', config.OUTPUT_PATH)
    
    # (Optional) Format graph if needed.
    # if config.ASSEMBLER == 'megahit':
    #     formatted_graph = graph_utils.format_graph(config.GRAPH_FILE_PATH)
    # elif config.ASSEMBLER == 'metaspades':
    #     formatted_graph = config.GRAPH_FILE_PATH  # or use formatted version
    # else:
    #     print('assembler not recognized, terminating...')
    #     return

    # # Create sequence file from graph.
    # seq_dict = graph_utils.create_seq_file_from_graph(formatted_graph, config.OUTPUT_PATH)
    
    # # Create graph from sequences.
    # g = graph_utils.create_graph(seq_dict, formatted_graph)
    # print(f"Graph built with {g.number_of_nodes()} nodes and {g.number_of_edges()} edges")
    
    # # Save graph if needed.
    # nx.write_gexf(g, os.path.join(config.OUTPUT_PATH, 'nx_graph.gexf'))
    g = nx.read_gexf(config.OUTPUT_PATH + 'nx_graph.gexf')
    print('graph read with ' + str(g.number_of_nodes()) + ' nodes and ' + str(g.number_of_edges()) + ' edges')

    # create_filtered_seq_file(g1)

    # Annotate ARGs.
    ARG_dict, ARG_target_dict = annotation.annotate_ARG()
    print("Total ARGs:", len(ARG_dict))
    with open(os.path.join(config.OUTPUT_PATH, 'ARG_dict_content.json'), 'w') as f:
        json.dump(ARG_dict, f, indent=4)
    
    ARG_dict = dict(sorted(ARG_dict.items(), key=lambda item: len(item[1])))

    # Extract ARG paths.
    ARG_path_dict, ARG_path_list = path_extraction.extract_ARG_paths_seq_gap(
        g, ARG_dict, config.DB_PATH_FASTA, config.OVERLAP, config.GENE_COV)
    
    # Optionally remove previous result files.
    for fname in ['left_seq.fasta', 'right_seq.fasta', 'contexts.txt', 'ARG_seq.fasta']:
        fpath = os.path.join(config.OUTPUT_PATH, f"result_contexts_cov_{config.GENE_COV}_radius_{config.RADIUS}", fname)
        if os.path.exists(fpath):
            os.remove(fpath)
    
    # Extract genomic contexts.
    context_extraction.extract_ARG_contexts(g, ARG_dict, ARG_path_dict, config.RADIUS, config.OVERLAP)
    print("Gene context extraction done.")

if __name__ == '__main__':
    main()
