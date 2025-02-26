# src/annotation.py
import os
import pandas as pd
from src.diamond_filter import filter_diamond_out

def annotate_ARG(db_path, output_path):
    """
    Runs diamond blastx against the deepARG-DB, filters the output, 
    and creates ARG_dict and ARG_target_dict from the filtered results.
    """
    # Run diamond blastx command
    # diamond_cmd = (
    #     f"diamond blastx -d {db_path} -q {output_path}generated_seq_file.fasta "
    #     f"--id 80 > {output_path}diamond_out"
    # )
    # os.system(diamond_cmd)
    
    # Filter the diamond output
    filter_diamond_out(f"{output_path}diamond_out", f"{output_path}filtered_diamond_out.tsv")  # This function should read diamond_out and write filtered_diamond_out.tsv

    # Read the filtered diamond output
    filtered_file = os.path.join(output_path, "filtered_diamond_out.tsv")
    df = pd.read_csv(filtered_file, delimiter='\t', header=None)
    df['Gene_Name'] = df[1].apply(lambda x: '|'.join(x.split('|')[1:]) if '|' in x else x)
    
    # Build ARG_dict and ARG_target_dict from the dataframe
    ARG_dict = df.groupby('Gene_Name').apply(lambda group: {
        node: max(
            [(min(start, end), max(start, end)) 
             for start, end in zip(group[group[0] == node][6], group[group[0] == node][7])],
            key=lambda se: se[1] - se[0]
        ) for node in group[0].unique()
    }).to_dict()
    
    ARG_target_dict = df.groupby('Gene_Name').apply(lambda group: {
        node: max(
            [(min(start, end), max(start, end)) 
             for start, end in zip(group[group[0] == node][8], group[group[0] == node][9])],
            key=lambda se: se[1] - se[0]
        ) for node in group[0].unique()
    }).to_dict()
    
    return ARG_dict, ARG_target_dict
