# src/path_extraction.py

from collections import defaultdict
from Bio import SeqIO
import networkx as nx 

def make_connected_groups(nodes, g1, total_gene_len, group_cov_ratio, max_node_dist):
    """
    makes groups with nodes(representing one ARG) which are connected with one another
    checks coverage for each group
    filters groups based on total coverage
    where MAX_NODE_DIST=10
    if one node doesn't have any other node from the nodes in its neighborhood,
    it belongs to the isolated_nodes
    returns updated nodes dict with only the nodes that belong to a group which covers the gene significantly
    """
    node_groups = {}
    group_assignments = {}
    group_id = 0
    isolated_nodes = []

    for node1 in list(nodes.keys()):
        is_isolated = True
        for node2 in list(nodes.keys()):
            if node1[:-1] != node2[:-1] and nx.has_path(g1, node1, node2):
                try:
                    dist = nx.shortest_path_length(g1, source=node1, target=node2)
                    if dist <= max_node_dist:
                        is_isolated = False
                        group1 = group_assignments.get(node1)
                        group2 = group_assignments.get(node2)

                        if group1 is None and group2 is None:
                            group_id += 1
                            node_groups[group_id] = {node1, node2}
                            group_assignments[node1] = group_id
                            group_assignments[node2] = group_id
                        elif group1 is not None and group2 is None:
                            node_groups[group1].add(node2)
                            group_assignments[node2] = group1
                        elif group1 is None and group2 is not None:
                            node_groups[group2].add(node1)
                            group_assignments[node1] = group2
                        elif group1 is not None and group2 is not None and group1 != group2:
                            node_groups[group1].update(node_groups[group2])
                            for node in node_groups[group2]:
                                group_assignments[node] = group1
                            del node_groups[group2]
                except nx.NetworkXNoPath:
                    pass
        
        if is_isolated:
            isolated_nodes.append(node1)
    
    for node in list(nodes.keys()):
        if node not in group_assignments:
            group_id += 1
            node_groups[group_id] = {node}
            group_assignments[node] = group_id

    filter_node_groups = {}
    for group_id, node_group in node_groups.items():
        covered_region = 0
        for node in node_group:
            covered_region += nodes[node][1] - nodes[node][0]
        if (covered_region / total_gene_len) >= group_cov_ratio:
            # node_groups.pop(group_id)
            filter_node_groups[group_id] = node_group
    print('filtered node groups ', filter_node_groups)

    node_list = set().union(*filter_node_groups.values())
    filtered_dict = {key: nodes[key] for key in node_list}

    return filtered_dict, isolated_nodes

def get_gene_path(ARG, nodes, g1, ref_ARGDB_dict, min_length_dict, overlap, max_gap_length, gene_cov, max_gap_node, group_cov_ratio, max_node_dist):

    """
    extracts the paths representing each ARG in the assembly graph and additional coverage info
    inputs: 
        ARG: specific gene
        nodes: nodes from the graph that aligns with that particular ARG
        g1: assembly graph in nx format
        ref_ARGDB_dict: dict built from reference ARG-DB
    outputs:
    filtered_unique_paths: list of paths representing the gene
    additional_path_info: ARG, path, covered_len: #bases aligned from the path to the gene, gene_len, covered_len/gene_len
    """
    print('current ARG ', ARG)
    print('number of nodes ', len(nodes))
    
    # gene_len = len(ref_ARGDB_dict[ARG].seq)*3
    gene_len = (min_length_dict[ARG])*3
    print('gene len ', gene_len)
    nodes, isolated_nodes = make_connected_groups(nodes, g1, gene_len, group_cov_ratio, max_node_dist)
    # isolated_nodes = find_isolated_nodes(g1, nodes)
    print('isolated nodes ', isolated_nodes)

    paths = []
    existing_node_set = set()
    covered_len_dict = {}

    def dfs(current_path, current_gap_len, covered_len, last_ARG_node, gap_node_count):

        # if ARG_count == len(node_groups[group_id]):
        #     paths.append(current_path)
        #     covered_len_dict[tuple(current_path)] = covered_len
        #     existing_node_set.update([node[:-1] for node in current_path])
        #     return

        last_node = current_path[-1]
        if last_node not in isolated_nodes:
            for neighbor in g1.neighbors(last_node):
                if neighbor not in current_path:
                    if neighbor in nodes: # neighbor is ARG
                    # if neighbor in node_groups[group_id]:
                        gap_len = current_gap_len + nodes[neighbor][0] - overlap
                        # if gap_len < MAX_GAP_LENGTH and target_nodes[last_ARG_node][1] <= target_nodes[neighbor][0]:
                        if gap_len < max_gap_length:
                            new_path = current_path + [neighbor]
                            covered_len += abs(nodes[neighbor][0] - nodes[neighbor][1])
                            dfs(new_path, (g1.nodes[neighbor]['length'] - nodes[neighbor][1]), covered_len, neighbor, 0)
                            
                    else: # neighbor is not ARG
                        gap_len = current_gap_len + g1.nodes[neighbor]['length'] - overlap
                        if gap_len < max_gap_length and gap_node_count < max_gap_node:
                            new_path = current_path + [neighbor]
                            dfs(new_path, gap_len, covered_len, last_ARG_node, gap_node_count + 1)
            
        if current_path[0] in nodes and current_path[-1] in nodes and not any(set(current_path) <= set(p) for p in paths):
            if covered_len / gene_len >= gene_cov:
                paths.append(current_path)
                covered_len_dict[tuple(current_path)] = covered_len
                existing_node_set.update([node[:-1] for node in current_path])

    # for group_id, node_group in node_groups.items():
    for node in nodes.keys():
        # for node in node_group:
        # for node, _ in nodes.items():
        if not node[:-1] in existing_node_set:
            initial_covered_len = abs(nodes[node][0] - nodes[node][1])
            dfs([node], (g1.nodes[node]['length'] - nodes[node][1]), initial_covered_len, node, 0)
        # else:
        #     paths.append(list(node_group))
        #     (node,) = node_group
        #     covered_len_dict[tuple(list(node_group))] = abs(nodes[node][0] - nodes[node][1])
        #     existing_node_set.update([node[:-1]])

    unique_paths = [p for p in paths if not any(set(p) < set(other_path) for other_path in paths if p != other_path)]
    filtered_unique_paths = []
    additional_path_info = []

    for path in unique_paths:
        filtered_unique_paths.append(path)
        additional_path_info.append((ARG, path, covered_len_dict[tuple(path)], gene_len, covered_len_dict[tuple(path)]/gene_len))

    return ARG, filtered_unique_paths, additional_path_info

def extract_ARG_paths_seq_gap(g1, ARG_dict, db_fasta, overlap, gene_cov, max_gap_length, max_gap_node, group_cov_ratio, max_node_dist):
    """
    Extracts ARG paths from the graph.
    """
    # Build min_length_dict from db_fasta
    min_length_dict = defaultdict(lambda: float('inf'))
    for record in SeqIO.parse(db_fasta, "fasta"):
        gene_name = '|'.join(record.description.split('|')[1:])
        sequence_length = len(record.seq)
        if sequence_length < min_length_dict[gene_name]:
            min_length_dict[gene_name] = sequence_length

    ref_ARGDB_dict = SeqIO.to_dict(SeqIO.parse(db_fasta, "fasta"))
    ARG_path_dict = {}
    ARG_path_list_2 = []

    for ARG, nodes in ARG_dict.items():
        # Call get_gene_path and update ARG_path_dict and ARG_path_list_2
        # Example (assuming get_gene_path is implemented):
        ARG, paths, additional_info = get_gene_path(ARG, nodes, g1, ref_ARGDB_dict, min_length_dict, overlap, max_gap_length, gene_cov, max_gap_node, group_cov_ratio, max_node_dist)
        ARG_path_dict[ARG] = paths
        ARG_path_list_2.extend(additional_info)
        # Write paths to file or process further if needed
    return ARG_path_dict, ARG_path_list_2
