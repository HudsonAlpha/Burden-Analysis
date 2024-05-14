#!/usr/bin/env python3

import pandas as pd
import sys

# ensure the correct number of arguments is provided
if len(sys.argv) != 4:
    print("Usage: ./script.py <wgsa_snp_anno>.tsv.gz <wgsa_indel_anno>.tsv.gz <merged_output_name>.tsv.gz")
    sys.exit(1)

# assign command-line arguments to variables
snp_file_path = sys.argv[1]
indel_file_path = sys.argv[2]
output_file_path = sys.argv[3]

# define the columns you want to keep and in the order you want them to appear
desired_columns = [
    '#CHROM', 'POS', 'REF', 'ALT', 
    'GenoCanyon_score', 'RegulomeDB_score', 
    'fathmmXF_score', 'Eigen_phred'
]

# map the desired columns to actual file columns (adjust as needed)
column_mapping = {
    'chr': '#CHROM',
    'pos': 'POS',
    'ref': 'REF',
    'alt': 'ALT',
    'GenoCanyon_score': 'GenoCanyon_score',
    'RegulomeDB_score': 'RegulomeDB_score',
    'fathmm-XF_score': 'fathmmXF_score',
    'Eigen-phred': 'Eigen_phred'
}

score_columns = ['GenoCanyon_score', 'fathmmXF_score', 'Eigen_phred']

# function to read, select, rename, and reorder columns
def process_file(file_path):
    # read the gzipped TSV file
    df = pd.read_csv(file_path, compression='gzip', sep='\t', usecols=column_mapping.keys(), low_memory=False)
        
    # rename the selected columns according to the mapping
    df = df.rename(columns=column_mapping)

    # reorder the DataFrame according to desired_columns
    df = df[desired_columns]
    
    df['#CHROM'] = 'chr' + df['#CHROM'].astype(str)
    
    return df

# function to process indels to collapse score to the best score in the set
def process_scores(score_field):
    if score_field == '.':
        return score_field
    
    # split the score_field, convert to numeric, handling non-numeric with 'coerce'
    scores = pd.to_numeric(score_field.split('|'), errors='coerce')
    
    # convert scores to a pandas Series to use the dropna method
    scores_series = pd.Series(scores)
    
    # use dropna on the pandas Series instead of the numpy array
    valid_scores = scores_series.dropna()
    
    # return the max of valid_scores if not empty, else return '.'
    return valid_scores.max() if not valid_scores.empty else '.'

# function to process indels to collapse score to the best score in the set for the regulomedb scores    
def process_regulomedb_score(score_field):
    # define the order of RegulomeDB scores from highest to lowest regulatory potential
    score_order = ['1a', '1b', '1c', '1d', '1e', '1f', '2a', '2b', '2c', '3a', '3b', '4', '5', '6', '7']
    
    # split the score field by '|' and keep only valid scores
    scores = score_field.split('|')
    valid_scores = [score for score in scores if score in score_order]
    
    # if no valid scores, return '.'
    if not valid_scores:
        return '.'
    
    # sort the valid scores by their regulatory potential, based on their position in score_order
    sorted_scores = sorted(valid_scores, key=lambda score: score_order.index(score))
    
    # return the score with the highest regulatory potential (lowest in terms of sorting)
    return sorted_scores[0]

# process each file
snp_df = process_file(snp_file_path)
indel_df = process_file(indel_file_path)

# process the indels in the scores
for col in score_columns:
    indel_df[col] = indel_df[col].astype(str).apply(process_scores)
    
indel_df['RegulomeDB_score'] = indel_df['RegulomeDB_score'].astype(str).apply(process_regulomedb_score)

# merge the DataFrames
merged_df = pd.concat([snp_df, indel_df], ignore_index=True)

merged_df['POS'] = pd.to_numeric(merged_df['POS'])

merged_df.sort_values(by=['#CHROM', 'POS'], inplace=True)

# write the merged DataFrame to a new gzipped TSV file
merged_df.to_csv(output_file_path, sep='\t', index=False)