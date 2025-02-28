# src/config.py

GRAPH_FILE_PATH = '/projects/ciwars/ARG_context/assemblies/CAMI_low/assembly_graph.fastg'
DB_PATH = '/projects/ciwars/ARG_context/database/filtered_deepARGDB.dmnd'
DB_PATH_FASTA = '/projects/ciwars/ARG_context/database/filtered_deepARGDB.fasta'
OUTPUT_PATH = '/projects/ciwars/ARG_context/ARGContextProfiler/results/'
READ1 = '/projects/ciwars/ARG_context/CAMI_synthetic_dataset/CAMI_low/CAMI_low_reads-1.fq'
READ2 = '/projects/ciwars/ARG_context/CAMI_synthetic_dataset/CAMI_low/CAMI_low_reads-2.fq'

OVERLAP = 55
RADIUS = 2000
GENE_COV = 0.6
MAX_GAP_LENGTH = 1000
MAX_NODE_DIST = 10
MAX_GAP_NODE = 10
MIN_CONTEXT_LEN = 100
GROUP_COV_RATIO = 0.4
ASSEMBLER = 'metaspades'
