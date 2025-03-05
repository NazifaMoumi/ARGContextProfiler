#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * This Nextflow pipeline implements the ARGContextProfiler pipeline.
 * It uses processes for each major step:
 *  - runAssembly: Runs metaspades to generate an assembly if no graph is provided.
 *  - useExistingGraph: Uses the provided graph if it exists.
 *  - buildDiamondDB: Builds the Diamond database.
 *  - createAssemblyGraph: Creates a NetworkX graph from the assembly graph.
 *  - runAnnotation: Runs the annotation step.
 *  - cleanOldResults: Cleans up previous result files.
 *  - runPathExtraction: Extracts ARG paths.
 *  - runContextExtraction: Extracts genomic contexts.
 *  - clusterContexts: Clusters contexts using mmseqs.
 *  - mapReads: Maps reads to the contexts using BWA/Samtools.
 *  - runReadAnalysis: Runs additional read analysis.
 *
 * Parameters are passed via Nextflow's params block.
 */

// params.graph       = params.graph       ?: '/projects/ciwars/ARG_context/assemblies/CAMI_low/assembly_graph.fastg'
// params.db_fasta    = params.db_fasta    ?: '/projects/ciwars/ARG_context/database/filtered_deepARGDB.fasta'
// params.output      = params.output      ?: '/projects/ciwars/ARG_context/ARGContextProfiler/CAMI_low_results/'
// params.read1       = params.read1       ?: '/projects/ciwars/ARG_context/CAMI_synthetic_dataset/CAMI_low/CAMI_low_reads-1.fq'
// params.read2       = params.read2       ?: '/projects/ciwars/ARG_context/CAMI_synthetic_dataset/CAMI_low/CAMI_low_reads-2.fq'
// params.overlap     = params.overlap     ?: 55
// params.radius      = params.radius      ?: 1000
// params.gene_cov    = params.gene_cov    ?: 0.6
// params.max_gap_length = params.max_gap_length ?: 1000
// params.max_node_dist   = params.max_node_dist   ?: 10
// params.max_gap_node    = params.max_gap_node    ?: 10
// params.min_context_len = params.min_context_len ?: 100
// params.group_cov_ratio = params.group_cov_ratio ?: 0.4
// params.assembler   = params.assembler   ?: 'metaspades'


params.graph       = '/projects/ciwars/ARG_context/assemblies/CAMI_low/assembly_graph.fastg'
params.db_fasta    = '/projects/ciwars/ARG_context/database/filtered_deepARGDB.fasta'
params.output      = 'CAMI_low_results'
params.read1       = '/projects/ciwars/ARG_context/CAMI_synthetic_dataset/CAMI_low/CAMI_low_reads-1.fq'
params.read2       = '/projects/ciwars/ARG_context/CAMI_synthetic_dataset/CAMI_low/CAMI_low_reads-2.fq'
params.overlap     = 55
params.radius      = 1000
params.gene_cov    = 0.6
params.max_gap_length = 1000
params.max_node_dist   = 10
params.max_gap_node    = 10
params.min_context_len = 100
params.group_cov_ratio = 0.4
params.assembler   = 'metaspades'

// // Channel for the input graph file (as a file object)
// graphInputCh = Channel.fromPath(params.graph, checkIfExists: false)

process runAssembly {
    tag "metaspades assembly"
    publishDir "${projectDir}/${output}/assemblies/", mode: 'copy'

    input:
        val output
        path read1
        path read2
    output:
        path "assemblies/assembly_graph.fastg", emit: assembly_graph
    script:
    """
    echo "Input graph not found. Running metaspades assembly..."
    mkdir -p ${projectDir}/${output}/assemblies
    metaspades.py -o assemblies/ -1 ${read1} -2 ${read2} --only-assembler -t 64 -m 400
    """
    // cp ${params.output}/assemblies/assembly_graph.fastg assembly_graph.fastg
}

process useExistingGraph {
    tag "existing graph"
    publishDir "${projectDir}/${output}", mode: 'copy'
    input:
        val output
    output:
        // file("assembly_graph.fastg") into graphCh
        path "assembly_graph.fastg", emit: assembly_graph
    script:
    """
    echo "Using existing graph file..."
    cp ${params.graph} assembly_graph.fastg
    """
}

process buildDiamondDB {
    tag "Diamond DB"
    publishDir "${projectDir}/${output}", mode: 'copy'

    input:
        path db_fasta
        val output
    output:
        path "query_db.dmnd", emit: query_db
    script:
    """
    echo "Building DIAMOND DB..."
    diamond makedb --in ${db_fasta} -d query_db
    """
}

process createAssemblyGraph {
    tag "Assembly Graph"
    publishDir "${projectDir}/${output}", mode: 'copy'

    input:
        path assembly_graph
        val output
        val assembler
    output:
        path "nx_graph.gexf", emit: nx_graph
        path "generated_seq_file.fasta", emit: seq_file
    script:
    """
    echo "Creating assembly graph..."
    python - <<EOF
from src import graph_utils
import os, networkx as nx
graph_file = "${assembly_graph}"
output = "${output}"
assembler = "${assembler}"
if assembler == 'megahit':
    formatted_graph = graph_utils.format_graph(graph_file)
elif assembler == 'metaspades':
    formatted_graph = graph_file
else:
    raise ValueError("Assembler not recognized: " + assembler)
seq_dict = graph_utils.create_seq_file_from_graph(formatted_graph, output)
g = graph_utils.create_graph(seq_dict, formatted_graph)
nx.write_gexf(g, "nx_graph.gexf")
EOF
    """
        // cp ${output}/nx_graph.gexf nx_graph.gexf
}

process runAnnotation {
    tag "Annotation"
    publishDir "${projectDir}/${output}", mode: 'copy'

    input:
        val output
        path query_db
        path seq_file
    output:
        path "ARG_dict_content.json", emit: arg_dict
        path "diamond_out", emit: diamond_out
        path "filtered_diamond_out.tsv", emit: filtered_diamond_out
    script:
    """
    echo "Annotating ARGs..."
    python - <<EOF
from src import annotation
import json, os
ARG_dict, ARG_target_dict = annotation.annotate_ARG("${query_db}", "${seq_file}")
with open(os.path.join("ARG_dict_content.json"), "w") as f:
    json.dump(ARG_dict, f, indent=4)
EOF
    """
}

process cleanOldResults {
    tag "Clean Old Results"
    input:
        file(arg_dict) from argDictCh
    script:
    """
    echo "Cleaning old result files..."
    rm -f ${params.output}/result_contexts_cov_${params.gene_cov}_radius_${params.radius}/left_seq.fasta
    rm -f ${params.output}/result_contexts_cov_${params.gene_cov}_radius_${params.radius}/right_seq.fasta
    rm -f ${params.output}/result_contexts_cov_${params.gene_cov}_radius_${params.radius}/contexts.txt
    rm -f ${params.output}/result_contexts_cov_${params.gene_cov}_radius_${params.radius}/ARG_seq.fasta
    """
}

process runPathExtraction {
    tag "Path Extraction"
    publishDir "${projectDir}/${output}", mode: 'copy'

    input:
        val output
        path nx_graph
        path arg_dict
    output:
        path "ARG_path_dict.json", emit: arg_path_dict
    script:
    """
    echo "Extracting ARG paths..."
    python - <<EOF
from src import path_extraction
import json, networkx as nx
g = nx.read_gexf("${nx_graph}")
with open("${arg_dict}") as f:
    ARG_dict = json.load(f)
ARG_path_dict, _ = path_extraction.extract_ARG_paths_seq_gap(g, ARG_dict, "${params.db_fasta}", float("${params.group_cov_ratio}"), int("${params.max_node_dist}"), int("${params.overlap}"), int("${params.max_gap_length}"), int("${params.max_gap_node}") ,float("${params.gene_cov}"))
with open("ARG_path_dict.json", "w") as out:
    json.dump(ARG_path_dict, out)
EOF
    """
}

process runContextExtraction {
    tag "Context Extraction"
    publishDir "${projectDir}/${output}", mode: 'copy'

    input:
        val output
        path nx_graph
        path arg_dict
        path arg_path_dict
    output:
        path "result_contexts_cov_${params.gene_cov}_radius_${params.radius}", emit: result_context_dir
    script:
    """
    echo "Extracting genomic contexts..."
    python - <<EOF
from src import context_extraction
import json, networkx as nx
g = nx.read_gexf("${nx_graph}")
with open("${arg_dict}") as f:
    ARG_dict = json.load(f)
with open("${arg_path_dict}") as f:
    ARG_path_dict = json.load(f)
context_extraction.extract_ARG_contexts(g, ARG_dict, ARG_path_dict, int("${params.radius}"), int("${params.overlap}"), int("${params.min_context_len}"), "${params.output}", float("${params.gene_cov}"))
EOF
    touch context_extraction.done
    """
}

process clusterContexts {
    tag "Cluster Contexts"
    publishDir "${projectDir}/${output}", mode: 'copy'

    input:
        val output
        path result_context_dir
    output:
        path "whole_context_with_length_clustered_rep_seq.fasta", emit: context_clustered
    script:
    """
    echo "Clustering contexts..."
    mmseqs easy-cluster ${result_context_dir}/whole_context_with_length.fasta \
      whole_context_with_length_clustered \
      tmp \
      --min-seq-id 0.95 -c 0.95 --cov-mode 1
    """
}

process mapReads {
    tag "Map Reads"
    publishDir "${projectDir}/${output}", mode: 'copy'

    input:
        val output
        path context_clustered
        path read1
        path read2
    output:
        path "read_mapping.bam", emit: read_mapping
    script:
    """
    echo "Mapping reads..."
    bwa index ${context_clustered}
    bwa mem -a -t 8 ${context_clustered} ${read1} ${read2} | samtools view -h -F 4 -b | samtools sort > read_mapping.bam
    samtools index read_mapping.bam
    """
}

process runReadAnalysis {
    tag "Read Analysis"
    publishDir "${projectDir}/${output}", mode: 'copy'

    input:
        val output
        path read_mapping
        path context_clustered
    output:

    script:
    """
    echo "Running read analysis..."
    python - <<EOF
from src import context_refinement
context_refinement.run_read_analysis("${read_mapping}", "${context_clustered}")
EOF
    """
}

workflow {

    def outDir = new File("${projectDir}/${params.output}")    
    if( outDir.exists() ) {
        println "Directory ${params.output} exists. Deleting its contents..."
        outDir.deleteDir()  // recursively delete the directory and its contents
    }
    outDir.mkdirs()

    // If input graph exists, use it; else run assembly
    if( !file(params.graph).exists() ) {
        assemblyGraphCh = runAssembly(params.output, params.read1, params.read2)
    } else {
        assemblyGraphCh = useExistingGraph(params.output)
    }

    diamondDBCh = buildDiamondDB(params.db_fasta, params.output)
    nxGraphCh = createAssemblyGraph(assemblyGraphCh.assembly_graph, params.output, params.assembler)
    argAnnotationCh = runAnnotation(params.output, diamondDBCh.query_db, nxGraphCh.seq_file)
    // cleanOldResults()
    pathExtractionCh = runPathExtraction(params.output, nxGraphCh.nx_graph, argAnnotationCh.arg_dict)
    contextExtractionCh = runContextExtraction(params.output, nxGraphCh.nx_graph, argAnnotationCh.arg_dict, pathExtractionCh.arg_path_dict)
    contextClusterCh = clusterContexts(params.output, contextExtractionCh.result_context_dir)
    readMapCh = mapReads(params.output, contextClusterCh.context_clustered, params.read1, params.read2)
    runReadAnalysis(params.output, readMapCh.read_mapping, contextClusterCh.context_clustered)

    println "Pipeline completed successfully! Output files located at ${params.output}"
}
