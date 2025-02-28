#!/usr/bin/env python

import os
import json
import networkx as nx
from src import config as default_config
from src import graph_utils, diamond_filter, annotation, path_extraction, context_extraction, context_refinement
import argparse
from dataclasses import dataclass
import subprocess
import shutil

@dataclass
class PipelineConfig:
    graph: str
    db_fasta: str
    output: str
    read1: str
    read2: str
    overlap: int
    radius: int
    gene_cov: float
    max_gap_length: int
    max_node_dist: int
    max_gap_node: int
    min_context_len: int
    group_cov_ratio: float
    assembler: str

def existing_path(path: str) -> str:
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(f"Path does not exist: {path}")
    return path

def parse_arguments() -> PipelineConfig:
    parser = argparse.ArgumentParser(description="Run ARG Context Pipeline")
    
    # File paths and directories
    parser.add_argument("--graph", default=default_config.GRAPH_FILE_PATH,
                        help="Path to the assembly graph in .fastg format")
    parser.add_argument("--db_fasta", type=existing_path, default=default_config.DB_PATH_FASTA,
                        help="Path to the reference DB FASTA file")
    parser.add_argument("--output", default=default_config.OUTPUT_PATH,
                        help="Directory for output results")
    parser.add_argument("--read1", type=existing_path, default=default_config.READ1,
                        help="Path to the read 1 of paired-end reads")
    parser.add_argument("--read2", type=existing_path, default=default_config.READ2,
                        help="Path to the read 2 of paired-end reads")

    # Parameters
    parser.add_argument("--overlap", type=int, default=default_config.OVERLAP,
                        help="Overlap length between nodes")
    parser.add_argument("--radius", type=int, default=default_config.RADIUS,
                        help="Radius for context extraction")
    parser.add_argument("--gene_cov", type=float, default=default_config.GENE_COV,
                        help="Minimum gene coverage ratio")
    parser.add_argument("--max_gap_length", type=int, default=default_config.MAX_GAP_LENGTH,
                        help="Maximum allowed gap length")
    parser.add_argument("--max_node_dist", type=int, default=default_config.MAX_NODE_DIST,
                        help="Maximum node distance")
    parser.add_argument("--max_gap_node", type=int, default=default_config.MAX_GAP_NODE,
                        help="Maximum allowed gap nodes")
    parser.add_argument("--min_context_len", type=int, default=default_config.MIN_CONTEXT_LEN,
                        help="Minimum context length")
    parser.add_argument("--group_cov_ratio", type=float, default=default_config.GROUP_COV_RATIO,
                        help="Group coverage ratio")
    parser.add_argument("--assembler", default=default_config.ASSEMBLER,
                        help="Assembler used (e.g., 'metaspades' or 'megahit')")
    
    args = parser.parse_args()
    return PipelineConfig(
        graph=args.graph,
        db_fasta=args.db_fasta,
        output=args.output,
        read1=args.read1,
        read2=args.read2,        
        overlap=args.overlap,
        radius=args.radius,
        gene_cov=args.gene_cov,
        max_gap_length=args.max_gap_length,
        max_node_dist=args.max_node_dist,
        max_gap_node=args.max_gap_node,
        min_context_len=args.min_context_len,
        group_cov_ratio=args.group_cov_ratio,
        assembler=args.assembler
    )

def run_assembly(cfg: PipelineConfig):
    """
    Checks if the input graph exists.
    If not, runs metaspades to generate the assembly graph.
    Updates cfg.graph to point to the new graph file.
    """
    if not os.path.exists(cfg.graph):
        print("Input graph not found. Running metaspades assembly to generate the graph...")
        assembly_dir = os.path.join(cfg.output, "assemblies")
        if not os.path.exists(assembly_dir):
            os.makedirs(assembly_dir)
        assembly_cmd = (
            f"metaspades.py -o {assembly_dir} -1 {cfg.read1} -2 {cfg.read2} "
            f"--only-assembler -t 64 -m 400"
        )
        subprocess.run(assembly_cmd, shell=True, check=True, stdout=subprocess.DEVNULL)
        new_graph = os.path.join(assembly_dir, "assembly_graph.fastg")
        if not os.path.exists(new_graph):
            raise FileNotFoundError("Assembly did not produce the expected graph file.")
        else:
            print(f"Assembly completed. New graph file located at {new_graph}")
            cfg.graph = new_graph

def build_diamond_db(cfg: PipelineConfig):
    """Create a DIAMOND database from the given FASTA file."""
    print("Building DIAMOND DB from the FASTA file...")
    db_path = os.path.join(cfg.output, "query_db")
    cmd = ["diamond", "makedb", "--in", cfg.db_fasta, "-d", db_path]
    subprocess.run(cmd, check=True)

def create_assembly_graph(cfg: PipelineConfig) -> nx.Graph:
    """Format the graph if needed, build the sequence file, and return the nx graph."""

    # Format if assembler is megahit
    if cfg.assembler == 'megahit':
        formatted_graph = graph_utils.format_graph(cfg.graph)
    elif cfg.assembler == 'metaspades':
        formatted_graph = cfg.graph
    else:
        raise ValueError(f"Assembler not recognized: {cfg.assembler}")

    # Create sequence file from graph
    seq_dict = graph_utils.create_seq_file_from_graph(formatted_graph, cfg.output)

    # Create and save the graph
    g = graph_utils.create_graph(seq_dict, formatted_graph)
    print(f"Graph built with {g.number_of_nodes()} nodes and {g.number_of_edges()} edges")
    nx.write_gexf(g, os.path.join(cfg.output, 'nx_graph.gexf'))

    # Reload the graph from GEXF (if you want to confirm or reuse)
    g = nx.read_gexf(os.path.join(cfg.output, 'nx_graph.gexf'))
    print(f"Graph read with {g.number_of_nodes()} nodes and {g.number_of_edges()} edges")
    return g

def run_annotation(cfg: PipelineConfig) -> dict:
    """Annotate ARGs using the pipeline's annotation module and DIAMOND."""
    ARG_dict, ARG_target_dict = annotation.annotate_ARG(cfg.output)
    print("Total ARGs:", len(ARG_dict))
    
    # Save to JSON
    with open(os.path.join(cfg.output, 'ARG_dict_content.json'), 'w') as f:
        json.dump(ARG_dict, f, indent=4)
    
    # Sort by number of nodes
    ARG_dict = dict(sorted(ARG_dict.items(), key=lambda item: len(item[1])))
    return ARG_dict

def run_path_extraction(g: nx.Graph, ARG_dict: dict, cfg: PipelineConfig) -> dict:
    """Extract ARG paths using path_extraction module."""
    ARG_path_dict, ARG_path_list = path_extraction.extract_ARG_paths_seq_gap(g, ARG_dict, cfg)
    return ARG_path_dict

def run_context_extraction(g: nx.Graph, ARG_dict: dict, ARG_path_dict: dict, cfg: PipelineConfig):
    """Extract genomic contexts."""
    context_extraction.extract_ARG_contexts(g, ARG_dict, ARG_path_dict, cfg)
    print("Gene context extraction done.")

def cluster_contexts(cfg: PipelineConfig):
    """Use mmseqs easy-cluster to reduce redundancy in the contexts."""
    print('Clustering the contexts to get rid of redundancy...')
    input_fasta = os.path.join(cfg.output, f"result_contexts_cov_{cfg.gene_cov}_radius_{cfg.radius}", "whole_context_with_length.fasta")
    out_prefix = os.path.join(cfg.output, f"result_contexts_cov_{cfg.gene_cov}_radius_{cfg.radius}", "whole_context_with_length_clustered")
    tmp_dir = os.path.join(cfg.output, f"result_contexts_cov_{cfg.gene_cov}_radius_{cfg.radius}", "tmp")
    
    cmd = [
        "mmseqs", "easy-cluster",
        input_fasta, out_prefix, tmp_dir,
        "--min-seq-id", "0.95", "-c", "0.95", "--cov-mode", "1"
    ]
    subprocess.run(cmd, check=True)

def map_reads(cfg: PipelineConfig):
    """Map input reads to the contexts for refinement using bwa + samtools."""
    print('Mapping input reads to the contexts for refinement...')
    rep_seq_fasta = os.path.join(cfg.output, f"result_contexts_cov_{cfg.gene_cov}_radius_{cfg.radius}", "whole_context_with_length_clustered_rep_seq.fasta")
    bam_out = os.path.join(cfg.output, f"result_contexts_cov_{cfg.gene_cov}_radius_{cfg.radius}", "read_mapping.bam")

    # bwa index
    subprocess.run(["bwa", "index", rep_seq_fasta], check=True, stdout=subprocess.DEVNULL)
    # bwa mem + samtools
    mem_cmd = [
        "bwa", "mem", "-a", "-t", "8",
        rep_seq_fasta, cfg.read1, cfg.read2
    ]

    full_cmd = (
        f"bwa mem -a -t 8 {rep_seq_fasta} {cfg.read1} {cfg.read2} "
        "| samtools view -h -F 4 -b "
        "| samtools sort "
        f"> {bam_out}"
    )
    subprocess.run(full_cmd, shell=True, check=True, stdout=subprocess.DEVNULL)

    # samtools index
    subprocess.run(["samtools", "index", bam_out], check=True, stdout=subprocess.DEVNULL)
    tmp_dir = os.path.join(cfg.output, f"result_contexts_cov_{cfg.gene_cov}_radius_{cfg.radius}", "tmp")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
        print("Temporary files removed.")

def clean_old_results(cfg: PipelineConfig):
    """Remove previous result files if they exist."""
    res_dir = os.path.join(cfg.output, f"result_contexts_cov_{cfg.gene_cov}_radius_{cfg.radius}")
    for fname in ["left_seq.fasta", "right_seq.fasta", "contexts.txt", "ARG_seq.fasta"]:
        fpath = os.path.join(res_dir, fname)
        if os.path.exists(fpath):
            os.remove(fpath)

def main():

    args = parse_arguments()

    print('Input graph:', args.graph)
    print('Output path:', args.output)

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    run_assembly(args)

    build_diamond_db(args)

    g = create_assembly_graph(args)

    ARG_dict = run_annotation(args)

    clean_old_results(args)

    # Extract ARG paths
    ARG_path_dict = run_path_extraction(g, ARG_dict, args)

    # Extract genomic contexts
    run_context_extraction(g, ARG_dict, ARG_path_dict, args)

    cluster_contexts(args)

    # Map reads for refinement
    map_reads(args)

    context_refinement.run_read_analysis(args)

    print(f"ARGContextProfiler completed successfully!\nOutput located at {args.output}")

if __name__ == '__main__':
    main()
