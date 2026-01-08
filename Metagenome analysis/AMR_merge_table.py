#!/usr/bin/env python3

"""
Requires python >= 3.9 and pandas and numpy
"""

import argparse
import os
import pandas as pd
import numpy as np 
import sys

parser = argparse.ArgumentParser(
                    description='Create TPM table of ARGs')
parser.add_argument('directory',type=str,
                    help='path to wgsa2 assembly outputs')

args = parser.parse_args()

merged_data = []

directory = args.directory

genes = 'genes/PREDgenes_ABUNtab.txt'

amr_info = 'genes/AMRs/AMRs_info.txt'

all_data=[]

final_data=[]

for sample in os.listdir(directory):

## import gene_TPM and AMR genes

	gene_TPM = pd.read_csv(os.path.join(directory, sample, genes), sep='\t')

	gene_TPM.to_csv("checking_TPM_file.csv")

	amr = pd.read_csv(os.path.join(directory, sample, amr_info), sep='\t')

	amrselect = amr[['Protein identifier', 'Accession of closest sequence']]

	amrselect.to_csv("select_check.tsv", sep="\t")


	amr_annotation = amr[['Gene symbol','Sequence name','Class','Subclass','Accession of closest sequence']]

	all_data.append(amr_annotation)

	amr_annotation_all = pd.concat(all_data, ignore_index=True)

	amr_annotation_all.to_csv("amr_annotations_all.tsv", sep="\t")
	
	amr_annotation_unique = amr_annotation_all.drop_duplicates()

	amr_annotation_na = amr_annotation_unique.fillna('Unknown')

	amr_annotation_na.to_csv("amr_annotations_unique.tsv", sep="\t")

	amr_annotation_agg = amr_annotation_na.groupby('Accession of closest sequence').agg({'Gene symbol': ', '.join, 'Sequence name': ', '.join,'Class': ', '.join,'Subclass': ', '.join}).reset_index()


	gene_TPMselect = gene_TPM[['#name','iTPM']]

	gene_TPMselect.to_csv("gene_TPMselect.tsv", sep="\t")

# merge data frames

	merged = amrselect.merge(gene_TPMselect, left_on='Protein identifier', right_on='#name')

	merged.to_csv("merged_amr_tpm_check.tsv", sep = "\t")

# remove unnecessary columns that might be making an issue with sum

	pruned = merged.drop(['Protein identifier', '#name'], axis=1)

	pruned.to_csv("check_removal.tsv", sep = "\t")

# move the Accession column to be the names of the rows 

	pruned.insert(0, 'Accession of closest sequence', pruned.pop('Accession of closest sequence'))

	pruned.to_csv("pruned_reorder_test.tsv", sep = "\t")

# sum by amr genes 
	sum = pruned.groupby(['Accession of closest sequence']).sum()

	sum_final = sum.rename(columns={'iTPM' : sample})

	final_data.append(sum_final)

merged_df = final_data[0]

for df in final_data[1:]:
	merged_df = pd.merge(merged_df, df, on='Accession of closest sequence', how='outer')

merged_df.fillna(0, inplace=True)

merged_df.to_csv("merged_data_raw.tsv", sep= "\t")

merged2 = pd.merge(merged_df, amr_annotation_agg, on='Accession of closest sequence', how="inner")


merged2.to_csv("Final_AMR_table.tsv", sep = "\t")


