#!/bin/python3
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import binom, poisson


exon_count_matrix_file = "all_rna_on_exon_reads.tsv"
transcriptome_fasta_file = "transcriptome.fa"
gene2transcript_file = "gene_id_transcript_id.tsv"

output_file = "pc_frequencies_mrl_poisson_fit.tsv"


replicates = ['rep1','rep2','rep3']

#sample name to fraction for MRL computation
sample_dict = {'CDS_RNAseq_Sample_HepG2polysomeFreeFraction' : 0,
               'CDS_RNAseq_Sample_HepG2polysome40s' : 0,
               'CDS_RNAseq_Sample_HepG2polysome80s' : 1,
               'CDS_RNAseq_Sample_HepG2polysomedisome' : 2,
               'CDS_RNAseq_Sample_HepG2polysometrisome' : 3,
               'CDS_RNAseq_Sample_HepG2polysometetrasome' : 4,
               'CDS_RNAseq_Sample_HepG2polysomepentasome' : 5,
               'CDS_RNAseq_Sample_HepG2polysomehexasome' : 6,
               'CDS_RNAseq_Sample_HepG2polysomeheptasome' : 7,
               'CDS_RNAseq_Sample_HepG2polysomeoctasome' : 8,
               'CDS_RNAseq_Sample_HepG2polysomenonasome' : 9,
               'CDS_RNAseq_Sample_HepG2polysomedecasome' : 10}
heavy_fract = 'CDS_RNAseq_Sample_HepG2polysomeheavysome'



'''
converts transcriptome fasta file to a pandas data frame 
with transcript ID and transcript length
returns data frame
takes name of the file as an input
'''
def extract_transcript_lengths_from_fasta(filename):
    len_dict = {}
    pc_dict = {}
    with open(filename) as f:
        transcript_id = "default"
        for line in f:
            if line[0] == ">":
                tmp = line.split(" ")
                transcript_id = tmp[0].strip().strip(">")
                len_dict[transcript_id] = 0

                if tmp[-1].strip().startswith("CDS="):
                    pc_dict[transcript_id] = True
                else:
                    pc_dict[transcript_id] = False
            else:
                len_dict[transcript_id] += len(line.strip())

    length_df = pd.DataFrame(len_dict.items(), columns=['transcript_ID', 'transcript_length'])
    pc_df = pd.DataFrame(pc_dict.items(), columns=['transcript_ID', 'protein_coding'])

    return pd.merge(length_df, pc_df, on='transcript_ID')



'''
computes transcript frequencies of genes in all fractions

requires an extra column with with transcript length called "transcript_length"
pools the counts over replicates

returns table with frequencies
'''
def compute_frequencies(df):
    df.set_index('gene_ID', inplace=True)

    #pool counts
    for fraction in [*sample_dict.keys(),heavy_fract]:
        df[fraction+'_pooled'] = df[[c for c in df.columns if fraction in c]].sum(axis=1)
    df = df[[*[c for c in df.columns if '_pooled' in c],'transcript_length']]

    normalized_df = df.div(df.transcript_length, axis=0)  #scale by length
    normalized_df.drop(columns=['transcript_length'], inplace=True)

    spikein_df = normalized_df[~normalized_df.index.str.startswith("ENSG")]  #split
    pc_df = normalized_df[normalized_df.index.str.startswith("ENSG")]

    spikein_norm_df = spikein_df.sum()  # compute normalization dfs
    pc_norm_df = pc_df.sum()

    frequencies = pc_df.div(spikein_norm_df.add(pc_norm_df, axis=0), axis=1)

    frequencies.reset_index(names=['gene_ID'], inplace=True)

    return frequencies



'''
define binomial distribution for fitting
'''
def binomial_distribution_function(x, n, p):
    return binom.pmf(x, n, p)



'''
define Poisson distribution for fitting
'''
def poisson_distribution_function(x, lam):
    return poisson.pmf(x, lam)



'''
determines chi^2 reduced
at:
	- fit function callable fit_func
	- fit parameter vector p
	- location vector x
	- function value vector y
	- error vector dy
'''
def chisqred(fit_func, p, x, y, dy):
    ret = 0.0

    for xi, yi, dyi in zip(x, y, dy):
        ret += ((fit_func(xi, *p) - yi) / dyi) ** 2

    dof = len(x) - len(p)

    if dof > 0:
        ret /= dof
    else:
        ret = np.nan

    return ret



'''
compute ribosome load
via fitting a binomial distribution
from input frequency series
assumes standard error of 0.01 for frequencies
'''
def compute_model_fit_RL(df):
    norm = df[heavy_fract+'_pooled']
    freq_dict = {}
    for name, fraction in sample_dict.items():
        tmp = df[name+'_pooled']

        norm += tmp
        if fraction in freq_dict:
            freq_dict[fraction] += tmp
        else:
            freq_dict[fraction] = tmp

    x_data = np.array(list(freq_dict.keys()))
    y_data = np.array(list(freq_dict.values())) / norm

    if not (x_data == np.arange(min(sample_dict.values()),max(sample_dict.values())+1)).all():
        print("Not all fractions specified are found in the input data")

    params, covmat = curve_fit(poisson_distribution_function, x_data, y_data, maxfev = 10000)

    return params[0], np.sqrt(np.diag(covmat))[0], chisqred(poisson_distribution_function,params,x_data,y_data,0.01*np.ones(len(x_data)))



'''
computes ribosome load (Sample et. al. Nature Biotech 2019)

input: data frame with frequencies of transcript in all fractions

returns data frame with ribosome load and error estimate
'''
def compute_ribosome_load(freqs):
    return freqs.apply(lambda x: compute_model_fit_RL(x), axis=1)




###############
#BEGIN SCRIPT
###############

# Load the dataset from the .tsv file, including the header
data = pd.read_csv(exon_count_matrix_file, sep="\t")
data = data.set_index('gene_ID')
data = data.reset_index(names=['gene_ID'])

# extract transcript lengths from fasta file
len_df = extract_transcript_lengths_from_fasta(transcriptome_fasta_file)

translate_df = pd.read_csv(gene2transcript_file, sep="\t", names=['gene_ID', 'transcript_ID'])
len_df = pd.merge(translate_df, len_df, on='transcript_ID')

len_df["transcript_length"] = len_df.groupby('gene_ID')["transcript_length"].transform('max')  #only keep longest
len_df = len_df.drop_duplicates(subset=["gene_ID"])

data = pd.merge(data, len_df[['gene_ID', 'transcript_length']], on="gene_ID")

frequencies = compute_frequencies(data)
frequencies = frequencies[np.count_nonzero(frequencies, axis=1) == len(frequencies.columns)]  #expression in all fractions

#compute MRL (first moment of Poissonian)
frequencies['ribosome_load'], frequencies['error_ribosome_load'], frequencies['relative_chisq_red'] = zip(*compute_ribosome_load(frequencies))

frequencies.to_csv("frequencies_mrl_poisson_fit.tsv",index=False,sep="\t")

#select protein coding
frequencies = pd.merge(frequencies, len_df[['gene_ID', 'protein_coding']], on='gene_ID')
frequencies = frequencies[frequencies['protein_coding']]
frequencies.drop(columns=["protein_coding"]).to_csv(output_file, index=False, sep="\t")

