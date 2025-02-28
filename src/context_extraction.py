# src/context_extraction.py

import os
from itertools import product

def extract_centered_paths(graph, center_node, fixed_length, ARG_path, overlap):
    centered_paths = []
    def dfs(node, current_path, current_length):
        nonlocal centered_paths
        if current_length == fixed_length or not list(graph.neighbors(node)):
            centered_paths.append(current_path)
            return
        for neighbor in graph.neighbors(node):
            if neighbor[:-1] in ARG_path:
                continue
            neighbor_length = graph.nodes[neighbor]['length']
            if neighbor not in current_path:
                if current_length + neighbor_length - overlap <= fixed_length:
                    dfs(neighbor, current_path + [neighbor], current_length + neighbor_length)
                else:
                    dfs(neighbor, current_path + [neighbor], fixed_length)
    for neighbor in list(graph.neighbors(center_node)):
        if neighbor[:-1] in ARG_path:
            continue
        neighbor_length = graph.nodes[neighbor]['length']
        if neighbor_length <= fixed_length:
            dfs(neighbor, [center_node, neighbor], neighbor_length)
        else:
            centered_paths.append([center_node, neighbor])
    return centered_paths

def write_results(ARG, ARG_seq, ARG_path, g1, left_seq_list, right_seq_list, all_paths, output_path, gene_cov, radius):
    
    if not os.path.exists(f"{output_path}result_contexts_cov_{gene_cov}_radius_{radius}/"):
        os.makedirs(f"{output_path}result_contexts_cov_{gene_cov}_radius_{radius}/")

    with open(f"{output_path}result_contexts_cov_{gene_cov}_radius_{radius}/left_seq.fasta", "a") as file:
        for idx, seq in enumerate(left_seq_list):
            header = f">{ARG}_{'_'.join(ARG_path)}_{idx}"
            file.write(header + "\n" + seq + "\n")

    with open(f"{output_path}result_contexts_cov_{gene_cov}_radius_{radius}/right_seq.fasta", "a") as file:
        for idx, seq in enumerate(right_seq_list):
            header = f">{ARG}_{'_'.join(ARG_path)}_{idx}"
            file.write(header + "\n" + seq + "\n")
    
    with open(f"{output_path}result_contexts_cov_{gene_cov}_radius_{radius}/contexts.txt", 'a') as file:
        for idx, path in enumerate(all_paths):
            header = f">{ARG}_{'_'.join(ARG_path)}_{idx}"
            file.write(header + "\n" + "->".join(path) + "\n")
    
    with open(f"{output_path}result_contexts_cov_{gene_cov}_radius_{radius}/ARG_seq.fasta", 'a') as file:
        header = f">{ARG}_{'_'.join(ARG_path)}"
        file.write(header + "\n" + ARG_seq + "\n")

    with open(f"{output_path}result_contexts_cov_{gene_cov}_radius_{radius}/whole_context_with_length.fasta", "a") as whole_file:
        for left_idx, left_seq in enumerate(left_seq_list):
            for right_idx, right_seq in enumerate(right_seq_list):
                header = f">{ARG}_{'_'.join(ARG_path)}_{left_idx}__{ARG}_{'_'.join(ARG_path)}_{right_idx}_length_"
                combined_seq = left_seq + ARG_seq + right_seq
                whole_file.write(header + str(len(combined_seq)) + "\n" + combined_seq + "\n")

def extract_merge_seq(all_paths, ARG, start_node, end_node, g1, ARG_dict, ARG_path, args):
    left_list = []
    right_list = []
    
    for path in all_paths:
        left_part = path[:path.index(start_node)]
        right_part = path[path.index(end_node) + 1:]
        
        if left_part:
            left_list.append(left_part)
        if right_part:
            right_list.append(right_part)

    output_list = []
    seen = set()

    for inner_list in left_list:
        tuple_inner = tuple(inner_list)
        if tuple_inner not in seen:
            seen.add(tuple_inner)
            output_list.append(inner_list)
    
    left_list = output_list
    
    output_list = []
    seen = set()

    for inner_list in right_list:
        tuple_inner = tuple(inner_list)
        if tuple_inner not in seen:
            seen.add(tuple_inner)
            output_list.append(inner_list)
    
    right_list = output_list
 
    start_pos, end_pos = ARG_dict[ARG][start_node][0], ARG_dict[ARG][end_node][1]
    rem_seq_right = g1.nodes[end_node]['sequence'][end_pos:]
    rem_seq_left = g1.nodes[start_node]['sequence'][:start_pos]
    
    right_seq_list = []
    left_seq_list = []
    
    if len(rem_seq_right) >= args.radius:
        rem_seq_right = rem_seq_right[0:args.radius]
        right_seq_list.append(rem_seq_right)
    else:
        if right_list:
            for right_nodes in right_list:
                right_seq = rem_seq_right
                for node in right_nodes:
                    seq2 = g1.nodes[node]['sequence']
                    # what if right_seq < overlap
                    right_seq = right_seq + seq2[args.overlap:]

                    if len(right_seq) >= args.radius:
                        right_seq = right_seq[0:args.radius]
                        break
                if len(right_seq) >= args.min_context_len:
                    right_seq_list.append(right_seq)
        else:
            if len(rem_seq_right) >= args.min_context_len:
                right_seq_list.append(rem_seq_right)
            
    if len(rem_seq_left) >= args.radius:
        rem_seq_left = rem_seq_left[-args.radius:]
        left_seq_list.append(rem_seq_left)
    else:
        if left_list:
            for left_nodes in left_list:
                left_seq = rem_seq_left

                for node in reversed(left_nodes):
                    seq2 = g1.nodes[node]['sequence']

                    left_seq = seq2[:-args.overlap] + left_seq

                    if len(left_seq) >= args.radius:
                        left_seq = left_seq[-args.radius:]
                        break
                if len(left_seq) >= args.min_context_len:
                    left_seq_list.append(left_seq)
        else:
            if len(rem_seq_left) >= args.min_context_len:
                left_seq_list.append(rem_seq_left)
    
    # ARG_path = all_paths[0][all_paths[0].index(start_node) : all_paths[0].index(end_node) + 1]
    
    if not ARG_path:
        print('testing')
        print(start_node)
        print(end_node)
        print(all_paths)

    ARG_seq = ''
    for node in ARG_path:
        if node == start_node and node == end_node:
            ARG_seq += g1.nodes[node]['sequence'][start_pos:end_pos]
        elif node == start_node:
            ARG_seq += g1.nodes[node]['sequence'][start_pos:]
        elif node == end_node:
            ARG_seq += g1.nodes[node]['sequence'][:end_pos]
        else:
            ARG_seq += g1.nodes[node]['sequence'][args.overlap:]
        
    write_results(ARG, ARG_seq, ARG_path, g1, left_seq_list, right_seq_list, all_paths, args.output, args.gene_cov, args.radius)


def extract_ARG_contexts(g1, ARG_dict, ARG_path_dict, args):
    mapping = {"+": "-", "-": "+"}
    for ARG, paths in ARG_path_dict.items():
        for path in paths:
            start_node, end_node = path[0], path[-1]
            start_pos, end_pos = ARG_dict[ARG][start_node][0], ARG_dict[ARG][end_node][1]
            if start_pos < args.radius:
                paths_minus = extract_centered_paths(g1, start_node[:-1] + mapping[start_node[-1]], args.radius, path, args.overlap)
            else:
                paths_minus = [[start_node[:-1] + mapping[start_node[-1]]]]
            if (g1.nodes[end_node]['length'] - end_pos) < args.radius:
                paths_plus = extract_centered_paths(g1, end_node, args.radius, path, args.overlap)
            else:
                paths_plus = [[end_node]]
            paths_minus = [list(reversed(pm)) for pm in paths_minus]
            paths_minus = [[contig[:-1] + mapping[contig[-1]] for contig in inner_lst] for inner_lst in paths_minus]
            all_paths = []
            if paths_minus and paths_plus:
                for pair in product(paths_minus, paths_plus):
                    merged_path = pair[0][:-1] + path + pair[1][1:]
                    all_paths.append(merged_path)
            elif paths_minus:
                for path_minus in paths_minus:
                    all_paths.append(path_minus[:-1] + path)
            elif paths_plus:
                for path_plus in paths_plus:
                    all_paths.append(path + path_plus[1:])
            else:
                print('No paths found for ', ARG, path)
                all_paths.append(path)
            extract_merge_seq(all_paths, ARG, start_node, end_node, g1, ARG_dict, path, args)
