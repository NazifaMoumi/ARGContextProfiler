# src/graph_utils.py

import os, re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import networkx as nx
import pyfastg

def create_seq_file_from_graph(formatted_graph_path, output_path):
    with open(formatted_graph_path, 'r') as graph_file:
        record_dict = SeqIO.to_dict(SeqIO.parse(graph_file, "fasta"),
                    key_function=lambda record: (record.id.split('_')[1] + ('-' if "'" in record.id.split(':')[0] else '+')))
    record_dict = {key: str(value.seq) for key, value in record_dict.items()}
    records = [SeqRecord(Seq(seq), id=seq_id, description="") for seq_id, seq in record_dict.items()]
    # out_file = os.path.join(output_path, 'generated_seq_file.fasta')
    out_file = os.path.join('generated_seq_file.fasta')
    with open(out_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
    return record_dict

def create_graph(seq_dict, formatted_graph_path):
    g = pyfastg.parse_fastg(formatted_graph_path)
    for u, v in g.edges():
        # Optionally, you can store edge lengths or other attributes if needed.
        pass
    nx.set_node_attributes(g, seq_dict, "sequence")
    return g

def delete_single_components(g):
    components = list(nx.weakly_connected_components(g))
    for comp in components:
        if len(comp) < 2:
            for node in comp:
                g.remove_node(node)
    return g

def create_filtered_seq_file(g, output_path):
    out_file = os.path.join(output_path, 'generated_seq_file.fasta')
    with open(out_file, 'w') as fw:
        for node in g.nodes():
            fw.write('>' + node + '\n' + g.nodes[node]['sequence'] + '\n')

def format_graph(graph_file_path):
    dir_path, filename = os.path.split(graph_file_path)
    base_filename, file_extension = os.path.splitext(filename)
    new_filename = f"{base_filename}_formatted{file_extension}"
    output_path = os.path.join(dir_path, new_filename)
    with open(graph_file_path, 'r') as infile, open(output_path, 'w') as outfile:
        records = SeqIO.parse(infile, "fasta")
        modified_records = []
        for record in records:
            header = record.id.replace("NODE", "EDGE")
            header = re.sub(r'_ID_\d+', '', header)
            record.id = header
            record.name = header
            record.description = ''
            modified_records.append(record)
        SeqIO.write(modified_records, outfile, "fasta")
    return output_path
