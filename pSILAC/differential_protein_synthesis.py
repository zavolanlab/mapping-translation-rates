#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import gc
import pysam
import json
import math
import statistics
import warnings
from scipy.optimize import curve_fit, OptimizeWarning


#file paths
intensity_filename = "All_Pep_Centric_peptide_centric_12012024.csv"

copy_numbers_filename = "HepG2_protein_copy_numbers.csv"
protein_IDs_filename = "trscrpt_ENSEMBL_gene_HGNC_uniprot.tsv"

time_cols_dict = { 0 : 'Log2Int_001_HP-0056_SLC_D23', \
	10 : 'Log2Int_004_HP-0056_SLC_D23', \
	30 : 'Log2Int_007_HP-0056_SLC_D23', \
	60 : 'Log2Int_010_HP-0056_SLC_D24', \
	120 : 'Log2Int_013_HP-0056_SLC_D23', \
	240 : 'Log2Int_019_HP-0056_SLC_D23', \
	360 : 'Log2Int_022_HP-0056_SLC_D23', \
	480 : 'Log2Int_025_HP-0056_SLC_D23', \
	660 : 'Log2Int_028_HP-0056_SLC_D23' }

replicate_ID = "rep1"

t_d = []
I_d = []


'''
intensity fit function
makes use of variable projection
to project the I_0 constant out
at point t with rate k^2
'''
def int_func(t, k):
	pSum = 0.0
	for i in range(len(t_d)):
		pSum += k*k*t_d[i] + I_d[i]
	return pSum / len(t_d) - k*k*t



'''
fit k and return it
'''
def fit_rate(raw_intensities):
	raw_times = np.fromiter(time_cols_dict.keys(),dtype=float)

	loc_t_d = []
	loc_I_d = []

	for i_t, inten in enumerate(raw_intensities):
		if inten > 0.0:
			loc_t_d.append(raw_times[i_t])
			loc_I_d.append(inten)

	global t_d
	t_d = loc_t_d
	global I_d
	I_d = loc_I_d	#I am deeply sorry to do this but fitting in python is incredibly pedestrian and this is probably the best way to achieve what I want

	if len(t_d) > 2:
		warnings.simplefilter("error",OptimizeWarning)
		try:
			popt, pcov = curve_fit(int_func, loc_t_d, loc_I_d, p0=np.sqrt(loc_I_d[0]/loc_I_d[-1]), maxfev=10000)
			return popt*popt * np.log(2.0)
		except OptimizeWarning:
			return -1.
	else:
		return -1.



'''
fits rates
for all peptides
'''
def fit_intensity_buildup_rates(df):
	rates = np.empty(shape=len(df))
	for index, row in df.iterrows():
		rates[index] = fit_rate(row[time_cols_dict.values()].to_numpy(dtype=float))
	ret_df = df.copy()
	ret_df['rate'] = rates
	ret_df = ret_df.dropna(subset=['rate'])
	ret_df = ret_df.sort_values(by=['protein_ID'])
	return ret_df



'''
average rates
over peptides
'''
def avg_rates(df):
	ret_df = pd.DataFrame(columns=['protein_ID','rate','spread'])
	ret_df['protein_ID'] = df['protein_ID'].unique()
	ret_df['rate'] = np.array(df.groupby('protein_ID')['rate'].mean())
	ret_df['spread'] = np.array(df.groupby('protein_ID')['rate'].max() - df.groupby('protein_ID')['rate'].min())
	return ret_df



'''
computes protein copy number per cell
averaged over three replicates
'''
def avg_prot_copy(df,columns):
	df[columns] = df[columns].replace(0, np.nan)
	df['rep_averaged_copy_number'] = df[columns].mean(axis=1)
	df[columns] = df[columns].replace(np.nan, 0)
	return df



'''
BEGIN script
'''
#load proteomics input file
print('Read input data.')
df_int = pd.read_csv(intensity_filename,usecols=['Protein',*time_cols_dict.values()])

#drop unmodified and ambiguous peptides
print('Drop unmodified peptides.')
df_int = df_int[(df_int['Protein'].str.contains("Lys8LFQ")) | (df_int['Protein'].str.contains("Arg10LFQ"))]

df_int = df_int[~(df_int['Protein'].str.contains(";"))]

#split the protein column into protein and peptide
df_int[['protein_ID','peptide']] = df_int['Protein'].str.split("_",expand=True)

df_int.reset_index(inplace=True)

#count_datapoints
df_int['nbr_data_points'] = np.count_nonzero(df_int[df_int.columns[df_int.columns.str.contains('Log2Int')]].to_numpy(),axis=1)

#do fitting
print('Do rate fitting per peptide.')
df_peptides_rates = fit_intensity_buildup_rates(df_int)

print('Done, print output to file.')
df_peptides_rates.to_csv('peptide_buildup_rates_'+replicate_ID+'.tsv',sep='\t',index=False)

df_peptides_rates = df_peptides_rates[df_peptides_rates['rate'] > 1e-5]
df_peptides_rates.to_csv('peptide_buildup_rates_tidy_'+replicate_ID+'.tsv',sep='\t',index=False)

print('Average over peptides')
df_protein_rates = avg_rates(df_peptides_rates)
df_protein_rates.to_csv('protein_buildup_rates_tidy_'+replicate_ID+'.tsv',sep='\t',index=False)


#load file with protein copy numbers
print('Normalize by protein abundance')
df_protein_copy_numbers = pd.read_csv(copy_numbers_filename,usecols=['Gene names','Averagecopy number H1','Averagecopy number H2','Averagecopy number H3'],decimal=",")

#remove duplicates
df_protein_copy_numbers = df_protein_copy_numbers.groupby('Gene names').mean().reset_index()

#translate to swissprot ids
translate_df = pd.read_csv(protein_IDs_filename,sep="\t")

#remove duplicates
translate_df = translate_df.drop_duplicates(subset=['swissprot_id'],keep='first',ignore_index=True)

#merge with translator df
df_protein_copy_numbers = pd.merge(df_protein_copy_numbers,translate_df,left_on="Gene names",right_on="HGNC_gene_symbol")
df_protein_copy_numbers = df_protein_copy_numbers[['swissprot_id','HGNC_gene_symbol','ensembl_gene_id','Averagecopy number H1','Averagecopy number H2','Averagecopy number H3']]

#compute average protein copy numbers
df_protein_copy_numbers = avg_prot_copy(df_protein_copy_numbers,['Averagecopy number H1','Averagecopy number H2','Averagecopy number H3'])

#merge with synthesis rates
df_protein_rates = pd.merge(df_protein_rates,df_protein_copy_numbers,left_on="protein_ID",right_on="swissprot_id")

#compute rates and spreads
df_protein_rates['protein_synthesis_rate'] = df_protein_rates['rate'] * df_protein_rates['rep_averaged_copy_number']
df_protein_rates['spread_protein_synthesis_rate'] = df_protein_rates['spread'] * df_protein_rates['rep_averaged_copy_number']

df_protein_rates = df_protein_rates[['protein_ID','HGNC_gene_symbol','ensembl_gene_id','rep_averaged_copy_number','protein_synthesis_rate','spread_protein_synthesis_rate']]
df_protein_rates = df_protein_rates[df_protein_rates['protein_synthesis_rate'] > 0.0]
df_protein_rates = df_protein_rates.fillna(0.0)

#save
df_protein_rates.to_csv('protein_synthesis_rates_'+replicate_ID+'.tsv',sep='\t',index=False)



