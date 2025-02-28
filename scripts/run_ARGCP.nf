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

params.graph       = params.graph       ?: '/projects/ciwars/ARG_context/assemblies/BB_2_2021_01_08_INF_S61_2/assembly_graph.fastg'
params.db_fasta    = params.db_fasta    ?: '/projects/ciwars/ARG_context/database/filtered_deepARGDB.fasta'
params.output      = params.output      ?: '/projects/ciwars/ARG_context/ARG-Context/hospital_2_analysis/'
params.read1       = params.read1       ?: '/projects/ciwars/ILLUMINA_decontaminated/Lane_3/BB_2_2021_01_08_INF_S61_R1_clean_decontam.fastq.gz'
params.read2       = params.read2       ?: '/projects/ciwars/ILLUMINA_decontaminated/Lane_3/BB_2_2021_01_08_INF_S61_R2_clean_decontam.fastq.gz'
params.overlap     = params.overlap     ?: 55
params.radius      = params.radius      ?: 2000
params.gene_cov    = params.gene_cov    ?: 0.6
params.max_gap_length = params.max_gap_length ?: 1000
params.max_node_dist   = params.max_node_dist   ?: 10
params.max_gap_node    = params.max_gap_node    ?: 10
params.min_context_len = params.min_context_len ?: 100
params.group_cov_ratio = params.group_cov_ratio ?: 0.4
params.assembler   = params.assembler   ?: 'metaspades'

// Channel for the input graph file (as a file object)
Channel.fromPath(params.graph, checkIfExists: false)
    .set { graphInputCh }

process runAssembly {
    tag "metaspades assembly"
    // This process runs only if the input graph doesn't exist.
    when:
        !file(params.graph).exists()
    output:
        file("assembly_graph.fastg") into graphCh
    script:
    """
    echo "Input graph not found. Running metaspades assembly..."
    mkdir -p ${params.output}/assemblies
    metaspades.py -o ${params.output}/assemblies/ -1 ${params.read1} -2 ${params.read2} --only-assembler -t 64 -m 400
    cp ${params.output}/assemblies/assembly_graph.fastg assembly_graph.fastg
    """
}

process useExistingGraph {
    tag "existing graph"
    when:
        file(params.graph).exists()
    output:
        file("assembly_graph.fastg") into graphCh
    script:
    """
    echo "Using existing graph file..."
    cp ${params.graph} assembly_graph.fastg
    """
}

process buildDiamondDB {
    tag "Diamond DB"
    input:
        file(db_fasta) from file(params.db_fasta)
    output:
        file("query_db.dmnd") into diamondDBCh
    script:
    """
    echo "Building DIAMOND DB..."
    diamond makedb --in ${db_fasta} -d ${params.output}/query_db
    cp ${params.output}/query_db.dmnd query_db.dmnd
    """
}

process createAssemblyGraph {
    tag "Assembly Graph"
    input:
        file(graph) from graphCh
    output:
        file("nx_graph.gexf") into graphGexfCh
    script:
    """
    echo "Creating assembly graph..."
    python - <<EOF
from src import graph_utils
import os, networkx as nx
graph_file = "${graph}"
output = "${params.output}"
seq_dict = graph_utils.create_seq_file_from_graph(graph_file, output)
g = graph_utils.create_graph(seq_dict, graph_file)
nx.write_gexf(g, os.path.join(output, "nx_graph.gexf"))
EOF
    """
}

process runAnnotation {
    tag "Annotation"
    input:
        file(graph_gexf) from graphGexfCh
    output:
        file("ARG_dict_content.json") into argDictCh
    script:
    """
    echo "Annotating ARGs..."
    python - <<EOF
from src import annotation
import json, os
output = "${params.output}"
ARG_dict, ARG_target_dict = annotation.annotate_ARG(output)
with open(os.path.join(output, "ARG_dict_content.json"), "w") as f:
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
    input:
        file(arg_dict) from argDictCh
        file(graph_gexf) from graphGexfCh
    output:
        file("ARG_path_dict.json") into argPathDictCh
    script:
    """
    echo "Extracting ARG paths..."
    python - <<EOF
from src import path_extraction
import json, networkx as nx
g = nx.read_gexf("${params.output}/nx_graph.gexf")
with open("${params.output}/ARG_dict_content.json") as f:
    ARG_dict = json.load(f)
ARG_path_dict, _ = path_extraction.extract_ARG_paths_seq_gap(g, ARG_dict, "${params.output}")
with open("ARG_path_dict.json", "w") as out:
    json.dump(ARG_path_dict, out)
EOF
    """
}

process runContextExtraction {
    tag "Context Extraction"
    input:
        file(arg_path_dict) from argPathDictCh
    output:
        file("context_extraction.done") into contextExtractionCh
    script:
    """
    echo "Extracting genomic contexts..."
    python - <<EOF
from src import context_extraction
import json, networkx as nx
g = nx.read_gexf("${params.output}/nx_graph.gexf")
with open("${params.output}/ARG_dict_content.json") as f:
    ARG_dict = json.load(f)
with open("ARG_path_dict.json") as f:
    ARG_path_dict = json.load(f)
context_extraction.extract_ARG_contexts(g, ARG_dict, ARG_path_dict, "${params.output}")
EOF
    touch context_extraction.done
    """
}

process clusterContexts {
    tag "Cluster Contexts"
    input:
        file(context_done) from contextExtractionCh
    output:
        file("whole_context_with_length_clustered_rep_seq.fasta") into clusteredFastaCh
    script:
    """
    echo "Clustering contexts..."
    mmseqs easy-cluster ${params.output}result_contexts_cov_${params.gene_cov}_radius_${params.radius}/whole_context_with_length.fasta \
      ${params.output}result_contexts_cov_${params.gene_cov}_radius_${params.radius}/whole_context_with_length_clustered \
      ${params.output}result_contexts_cov_${params.gene_cov}_radius_${params.radius}/tmp \
      --min-seq-id 0.95 -c 0.95 --cov-mode 1
    cp ${params.output}result_contexts_cov_${params.gene_cov}_radius_${params.radius}/whole_context_with_length_clustered_rep_seq.fasta .
    """
}

process mapReads {
    tag "Map Reads"
    input:
        file(clustered_fasta) from clusteredFastaCh
    output:
        file("read_mapping.bam") into readMappingBamCh
    script:
    """
    echo "Mapping reads..."
    bwa index ${clustered_fasta}
    bwa mem -a -t 8 ${clustered_fasta} ${params.read1} ${params.read2} | samtools view -h -F 4 -b | samtools sort > read_mapping.bam
    samtools index read_mapping.bam
    """
}

process runReadAnalysis {
    tag "Read Analysis"
    input:
        file(bam) from readMappingBamCh
    output:
        file("read_analysis.done") into readAnalysisDoneCh
    script:
    """
    echo "Running read analysis..."
    python - <<EOF
from src import context_refinement
context_refinement.run_read_analysis("${params.output}")
EOF
    touch read_analysis.done
    """
}

workflow {
    // If input graph exists, use it; else run assembly
    if( !file(params.graph).exists() ) {
        runAssembly()
    } else {
        useExistingGraph()
    }
    buildDiamondDB()
    createAssemblyGraph()
    runAnnotation()
    cleanOldResults()
    runPathExtraction()
    runContextExtraction()
    clusterContexts()
    mapReads()
    runReadAnalysis()

    println "Pipeline completed successfully! Output located at ${params.output}"
}
