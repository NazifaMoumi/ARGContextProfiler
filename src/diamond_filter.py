# src/diamond_filter.py

import csv
import pandas as pd

def filter_diamond_out(input_file, output_file):
    """
    Filters the diamond output to keep alignments with alignment_length >= 50 
    and the lowest e-value per group.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        for row in reader:
            if len(row) > 2 and int(row[3]) >= 50:
                writer.writerow(row)

    df = pd.read_csv(output_file, sep='\t', header=None,
                     names=['query_id', 'target_id', 'percentage_identity', 'alignment_length',
                            'mismatches', 'gap_opens', 'query_start', 'query_end',
                            'target_start', 'target_end', 'e_value', 'bit_score'])
    df['comparison_key'] = df.apply(lambda row: f"{row['query_id']}_{row['query_start'] // 100}_{row['query_end'] // 100}", axis=1)
    filtered_df = df.loc[df.groupby('comparison_key')['e_value'].idxmin()]
    filtered_df.drop(columns=['comparison_key'], inplace=True)
    filtered_df.to_csv(output_file, sep='\t', index=False, header=False)
