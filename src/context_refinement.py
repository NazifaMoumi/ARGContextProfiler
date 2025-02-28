# src/context_refinement.py

import os
import math
import pysam
import pandas as pd
import numpy as np
from Bio import SeqIO
import collections

# You can define or import these from config, if needed:
BOUNDARY = 500

def calculate_median_dev(cfg):
    samfile = pysam.AlignmentFile(cfg + 'read_mapping.bam', "rb")
    
    all_reads = samfile.fetch()
    size_freq = collections.defaultdict(int)
    for read in all_reads:
        if read.rnext == read.tid and read.is_paired:
            size = abs(read.isize)
            size_freq[size] += 1
            
    all_size = []
    for key, value in size_freq.items():
        all_size.extend([key] * int(value))
    median_size = np.median(all_size)
    residuals = abs(np.array(all_size) - median_size)
    mad_size = 1.4826 * np.median(residuals)
    return median_size, mad_size

def read_pair_consistency(cfg):
    """
    Analyzes read pair consistency using 'read_mapping.bam' in cfg.output.
    Writes 'read_pair_consistency_scores.csv' to cfg.output.
    """
    bam_path = os.path.join(cfg, 'read_mapping.bam')
    samfile = pysam.AlignmentFile(bam_path, 'rb')

    median_size, mad_size = calculate_median_dev(cfg)
    references = samfile.references
    lengths = samfile.lengths

    read_dicts = {
        "contig": [], "start_pos": [], "read_count": [], "proper_read_count": [],
        "inversion_read_count": [], "clipped_read_count": [],
        "supplementary_read_count": [], "discordant_size_count": [],
        "discordant_loc_count": [], "length": []
    }

    count_empty = 0
    excluded = 0

    for ref, lens in zip(references, lengths):
        contig_reads = list(samfile.fetch(ref.strip()))
        read_dict = {
            "start_pos": [], "read_count": [], "proper_read_count": [],
            "inversion_read_count": [], "clipped_read_count": [],
            "supplementary_read_count": [], "discordant_size_count": [],
            "discordant_loc_count": []
        }
        read_temp = {
            "num_read": 0, "num_proper": 0, "num_inversion": 0, "num_clipped": 0,
            "num_supplementary": 0, "num_discordant_size": 0, "num_discordant_loc": 0
        }
        pos = 0

        if len(contig_reads) == 0:
            count_empty += 1

        newpos_list = []
        for read in contig_reads:
            new_pos = math.floor((read.reference_start - BOUNDARY) / 100) * 100 + BOUNDARY
            newpos_list.append(new_pos)
            if read.reference_start < BOUNDARY:
                continue
            if pos == 0:
                pos = new_pos
            elif new_pos != pos:
                # Store the current bin's counts
                read_dict["start_pos"].append(pos)
                read_dict["read_count"].append(read_temp["num_read"])
                read_dict["proper_read_count"].append(read_temp["num_proper"])
                read_dict["inversion_read_count"].append(read_temp["num_inversion"])
                read_dict["clipped_read_count"].append(read_temp["num_clipped"])
                read_dict["supplementary_read_count"].append(read_temp["num_supplementary"])
                read_dict["discordant_size_count"].append(read_temp["num_discordant_size"])
                read_dict["discordant_loc_count"].append(read_temp["num_discordant_loc"])

                # Reset for new bin
                read_temp = {
                    "num_read": 0, "num_proper": 0, "num_inversion": 0, "num_clipped": 0,
                    "num_supplementary": 0, "num_discordant_size": 0, "num_discordant_loc": 0
                }
                pos = new_pos

            read_temp["num_read"] += 1

            if read.is_paired:
                if read.rnext == read.tid:  # same contig
                    if read.is_proper_pair:
                        read_temp["num_proper"] += 1
                    # If both reads in the pair have the same orientation
                    if (read.is_reverse + read.mate_is_reverse) != 1:
                        read_temp["num_inversion"] += 1
                    # Check if insert size is within median +/- 3*MAD
                    if not (median_size - 3 * mad_size <= abs(read.isize) <= median_size + 3 * mad_size):
                        read_temp["num_discordant_size"] += 1
                else:
                    # read pair is on a different contig
                    read_temp["num_discordant_loc"] += 1

            # Check for clipping
            if read.get_cigar_stats()[0][4] > 20:
                read_temp["num_clipped"] += 1
            # Check for supplementary
            if read.is_supplementary and read.get_cigar_stats()[0][5] > 20:
                read_temp["num_supplementary"] += 1

        # If all reads fall into a single bin
        if len(set(newpos_list)) == 1 and all(len(v) == 0 for v in read_dict.values()):
            excluded += 1
            read_dict["start_pos"].append(pos)
            read_dict["read_count"].append(read_temp["num_read"])
            read_dict["proper_read_count"].append(read_temp["num_proper"])
            read_dict["inversion_read_count"].append(read_temp["num_inversion"])
            read_dict["clipped_read_count"].append(read_temp["num_clipped"])
            read_dict["supplementary_read_count"].append(read_temp["num_supplementary"])
            read_dict["discordant_size_count"].append(read_temp["num_discordant_size"])
            read_dict["discordant_loc_count"].append(read_temp["num_discordant_loc"])
            read_temp = {
                "num_read": 0, "num_proper": 0, "num_inversion": 0, "num_clipped": 0,
                "num_supplementary": 0, "num_discordant_size": 0, "num_discordant_loc": 0
            }

        # Extend the master dict
        read_dicts["start_pos"].extend(read_dict["start_pos"])
        read_dicts["contig"].extend([ref] * len(read_dict["start_pos"]))
        read_dicts["read_count"].extend(read_dict["read_count"])
        read_dicts["proper_read_count"].extend(read_dict["proper_read_count"])
        read_dicts["inversion_read_count"].extend(read_dict["inversion_read_count"])
        read_dicts["clipped_read_count"].extend(read_dict["clipped_read_count"])
        read_dicts["supplementary_read_count"].extend(read_dict["supplementary_read_count"])
        read_dicts["discordant_size_count"].extend(read_dict["discordant_size_count"])
        read_dicts["discordant_loc_count"].extend(read_dict["discordant_loc_count"])
        read_dicts["length"].extend([lens] * len(read_dict["start_pos"]))

    print('Number of empty contexts:', count_empty)
    print('Excluded bins:', excluded)

    data = pd.DataFrame(read_dicts)
    # Drop the 'start_pos' column before grouping
    data.drop(columns=['start_pos'], inplace=True)

    # Group by contig and take mean
    data = data.groupby('contig', as_index=False).mean()

    # Split into left and right contigs
    data[['left', 'right']] = data['contig'].str.split('__', expand=True)
    data.drop(columns=['contig'], inplace=True)

    # Adjust the 'right' column to remove the last two tokens
    data['right'] = data['right'].str.split('_').str[:-2].str.join('_')

    # group_id is everything except the last token of 'left'
    data['group_id'] = data['left'].str.split('_').str[:-1].str.join('_')

    # Reorder columns
    data = data[[
        'group_id', 'left', 'right',
        'read_count', 'proper_read_count', 'inversion_read_count',
        'clipped_read_count', 'supplementary_read_count',
        'discordant_size_count', 'discordant_loc_count', 'length'
    ]]

    out_csv = os.path.join(cfg, 'read_pair_consistency_scores.csv')
    data.to_csv(out_csv, index=False)
    return data

def read_coverage_uniformity(cfg):
    """
    Analyzes coverage uniformity using 'read_mapping.bam' in cfg.output.
    Writes 'read_coverage_uniformity_scores.csv' to cfg.output.
    """
    bam_path = os.path.join(cfg, 'read_mapping.bam')
    samfile = pysam.AlignmentFile(bam_path, 'rb')
    
    references = samfile.references
    lengths = samfile.lengths

    frag_dict = {
        "contig": [], "start_pos": [],
        "normalized_fragment_coverage": [], "normalized_fragment_deviation": []
    }
    count_empty = 0

    for ref, lens in zip(references, lengths):
        reads = list(samfile.fetch(ref))
        if len(reads) == 0:
            count_empty += 1
        
        frag_coverage = np.zeros(lens, dtype=int)

        # Mark coverage for each read (or read pair)
        for read in reads:
            size = abs(read.isize)
            start = min(read.next_reference_start, read.reference_start, read.reference_end)
            if size > 0:
                end = start + size
            else:
                end = max(read.next_reference_start, read.reference_start, read.reference_end)
            if start < 0: 
                start = 0
            if end > lens:
                end = lens
            if start == end:
                continue
            frag_coverage[start:end] += 1
        
        cov = {"pos": [], "coverage": [], "deviation": []}
        # Slide through in bins of 100 ignoring boundary
        for i in range(BOUNDARY, lens, 100):
            start_bin = i
            end_bin = i + 100
            if end_bin >= lens - BOUNDARY:
                break
            bin_cov = frag_coverage[start_bin:end_bin]
            mean_cov = np.mean(bin_cov) if len(bin_cov) > 0 else 0
            std_cov = np.sqrt(np.var(bin_cov)) if len(bin_cov) > 0 else 0
            deviation = std_cov / mean_cov if mean_cov != 0 else 0

            cov["pos"].append(start_bin)
            cov["coverage"].append(mean_cov)
            cov["deviation"].append(deviation)

        for i in range(len(cov["pos"])):
            frag_dict["contig"].append(ref)
            frag_dict["start_pos"].append(cov["pos"][i])
            frag_dict["normalized_fragment_coverage"].append(cov["coverage"][i])
            frag_dict["normalized_fragment_deviation"].append(cov["deviation"][i])

    print('Number of empty contexts (no reads):', count_empty)
    data = pd.DataFrame(frag_dict).fillna(0)
    
    # out_csv = os.path.join(cfg, 'read_coverage_uniformity.csv')
    # data.to_csv(out_csv, index=False)

    # # Now group by contig, compute mean coverage and deviation
    # data = pd.read_csv(out_csv)
    data.drop(columns=['start_pos'], inplace=True)
    data = data.groupby('contig', as_index=False)[
        ['normalized_fragment_coverage', 'normalized_fragment_deviation']
    ].mean()

    # Split into left and right contigs
    data[['left', 'right']] = data['contig'].str.split('__', expand=True)
    data.drop(columns=['contig'], inplace=True)

    # Adjust the 'right' column
    data['right'] = data['right'].str.split('_').str[:-2].str.join('_')

    # group_id is everything except the last token of 'left'
    data['group_id'] = data['left'].str.split('_').str[:-1].str.join('_')

    data = data[[
        'group_id', 'left', 'right',
        'normalized_fragment_coverage', 'normalized_fragment_deviation'
    ]]

    out_csv2 = os.path.join(cfg, 'read_coverage_uniformity_scores.csv')
    data.to_csv(out_csv2, index=False)
    return data

def filter_contexts_read_outliers(cfg):
    """
    Merges read_pair_consistency_scores.csv and read_coverage_uniformity_scores.csv,
    identifies outliers, and filters them out from 
    'whole_context_with_length_clustered_rep_seq.fasta'.
    """
    read_pair_file = os.path.join(cfg, 'read_pair_consistency_scores.csv')
    coverage_file = os.path.join(cfg, 'read_coverage_uniformity_scores.csv')
    merged_out = os.path.join(cfg, 'all_read_scores.csv')

    df_read_pair = pd.read_csv(read_pair_file)
    df_read_coverage = pd.read_csv(coverage_file)

    merged_df = pd.merge(df_read_pair, df_read_coverage, on=['group_id','left','right'], how='outer')
    merged_df.to_csv(merged_out, index=False)

    df = pd.read_csv(merged_out)
    features = [
        'read_count','proper_read_count','inversion_read_count','clipped_read_count',
        'supplementary_read_count','discordant_size_count','discordant_loc_count',
        'normalized_fragment_coverage','normalized_fragment_deviation'
    ]

    # Define which features are "higher is worse" (up_features) vs. "lower is worse" (dn_features)
    up_features = [
        'normalized_fragment_deviation','discordant_loc_count','discordant_size_count',
        'supplementary_read_count','clipped_read_count','inversion_read_count'
    ]
    dn_features = [
        'read_count','proper_read_count','normalized_fragment_coverage'
    ]

    def find_outliers(group, feature_list):
        is_outlier = pd.Series([False] * len(group), index=group.index)
        for feature in feature_list:
            mean_val = group[feature].mean()
            std_val = group[feature].std()
            if pd.isna(std_val):
                continue
            if feature in up_features:
                upper_bound = mean_val + 3*std_val
                is_outlier |= group[feature] > upper_bound
            elif feature in dn_features:
                lower_bound = mean_val - 3*std_val
                is_outlier |= group[feature] < lower_bound
        return group[is_outlier]

    all_feature_outliers_df = df.groupby('group_id').apply(lambda g: find_outliers(g, features))
    if not all_feature_outliers_df.empty:
        print("Outlier rows:\n", all_feature_outliers_df)
        print("Number of outlier rows:", len(all_feature_outliers_df))
        all_feature_outliers_df['combined'] = (
            all_feature_outliers_df['left'].astype(str) + '__' + all_feature_outliers_df['right'].astype(str)
        )
        remove_substrings = set(all_feature_outliers_df['combined'].tolist())
    else:
        print("No outliers found.")
        remove_substrings = set()

    # Filter out outliers from the FASTA
    rep_fasta = os.path.join(cfg, "whole_context_with_length_clustered_rep_seq.fasta")
    filtered_fasta = os.path.join(cfg, "whole_context_with_length_clustered_filtered.fasta")
    
    if not os.path.exists(rep_fasta):
        print(f"Warning: {rep_fasta} does not exist. Cannot filter contexts.")
        return

    filtered_sequences = []
    for record in SeqIO.parse(rep_fasta, 'fasta'):
        header = record.description
        if not any(sub in header for sub in remove_substrings):
            filtered_sequences.append(record)

    with open(filtered_fasta, 'w') as output_handle:
        SeqIO.write(filtered_sequences, output_handle, 'fasta')

    print("Filtered contexts written to:", filtered_fasta)

def run_read_analysis(cfg):
    """
    Runs all three read analysis steps in sequence:
      1. read_pair_consistency
      2. read_coverage_uniformity
      3. filter_contexts_read_outliers
    """
    read_mapping_path = f"{cfg.output}result_contexts_cov_{cfg.gene_cov}_radius_{cfg.radius}/"
    print("Running read pair consistency analysis...")
    read_pair_consistency(read_mapping_path)

    print("Analyzing coverage uniformity...")
    read_coverage_uniformity(read_mapping_path)

    print("Filtering outlier contexts based on read scores...")
    filter_contexts_read_outliers(read_mapping_path)

    print("Read analysis steps completed.")
