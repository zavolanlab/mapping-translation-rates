# Mapping translation rates
This is a repository containing the analysis scripts for our comprehensive mapping of translation rates in the HepG2 liver cancer cell line.
The underlying data set is a multi-omics data set with data from RNA sequencing, ribosome profiling, and proteomics.
The data encompasses data from widely used steady state techniques like RNA sequencing, ribosome profiling, and proteomics.
Moreover, we use RNA sequencing of different fractions obtained from polysome profiling to calculate mean ribosome load measures, ribosome profiling time courses after harringtonine treatment to observe ribosome runoff, and pulsed SILAC to measure the rate of protein buildup in HepG2 cells, i.e. we complement the standard (steady-state) techniques with methods that can resolve the time-dimension. While data analysis of the former is quite well-established, we publish the analysis scripts for the less standard methods (polysome profiling, ribosome runoff, pSILAC).

## Mean ribosome load from polysome profiling
The script ```polysome_analysis_model_RL_pool.py``` performs the analysis of polysome profiling data.

### Input data
The script needs three input files and some column names:
1. a count matrix of RNA seq reads mapped to any exon of a particular gene (tsv format)
2. a tsv file that coordinates gene IDs with transcript IDs (needs columns named "gene_ID" and "transcript_ID")
3. a fasta file of the transcriptome (possibly assemble with gffread) to extract the transcript length and if the transcript codes for a protein
4. specify column names of the fractions together with the corresponding number of ribosomes in the ```sample_dict``` dictionary, and the replicate names in the ```replicates``` array. The remainder fraction is treated separately and its name should be specified as the ```heavy_fract```.

The names of these input files can be specified in the upper paragraph of the script.

### Output data
The produced output is a tsv file containing the frequency of one gene in each fraction, the ribosome load, the fitting error of the ribosome load, and a *relative* estimate of the reduced chi^2. Relative means its absolute size is not meaningful, but it can be used as a criterion for fit likelihood when comparing genes.

### Procedure
Experimentally, an assessment of the total amount of RNA in each fraction is possible by processing the isolated RNA together with a given amount of spikein RNA, that is kept constant among fractions. To this end, the spikein sequences have to be amended to the genome that is used for mapping (ENSEMBL GRCH 38.105 and lexogen SIRV set 3 in our case).

Computationally, we extract mappings to exons with Rsubread for our initial count matrix. Then, the script:
- extracts transcript lengths and if they are coding for proteins from the transcriptome fasta file
- computes the frequencies of transcripts (per fraction, ratio of TPMs) from the count matrix
- computes relative representation of each fraction (per gene) and fits a Poisson distribution (first moment of the Poisson distribution is lambda, *i.e.* the mean ribosome load)
- saves the frequencies, lambda, error of lambda, and relative reduced chi^2 into a tsv file

This method is superior to just computing the weighted average of ribosomes per gene, since it is not intrinsically limited to the maximum fraction that you can resolve (10 in our case).

## Protein synthesis rates from pSILAC
The script ```differential_protein_synthesis.py``` performs the analysis of pSILAC data.

### Input data
The script needs the following inputs:
1. a .csv file with (log2-fold) peptide intensities
2. specify the columns in the ```time_cols_dict``` dictionary, together with the time differences (in minutes) between switching from heavy to light amino acids and measurement
3. the identifier of the replicate to analyze
4. a .csv file with the protein copy numbers per cell (e.g. Wisniewski *et. al.*. 2016. “In-Depth Quantitative Analysis and Comparison of the Human Hepatocyte and Hepatoma Cell Line HepG2 Proteomes.” Journal of Proteomics)
5. specify the column names over which the protein copy number should be averaged in ```copy_number_columns```
6. a .tsv file that connects ENSEMBL gene IDs to transcript IDs, swissprot and uniprot IDs

### Procedure
With this input at hand, we proceed as follows:
1. only keep the peptides that are labeled as heavy (in our case either by ```Lys8LFQ```or ```Arg10LFQ```)
2. read off the protein id and the peptide sequence from the peptide identifier
3. count the number of time points which are resolved per peptide
4. fit a linear curve with negative slope to the intensities of the heavily labeled peptides (positive slope if doing the transition light -> heavy at t=0). Since Peptide intensities themselves are not meaningful, we can project them out using variable projection (O’Leary, Rust 2012: Variable Projection for Nonlinear Least Squares Problems) to make the fit more stable.
5. average over multiple peptides being a part of a single protein to end up with protein buildup rates
6. normalize by protein abundance per cell to end up with protein synthesis rates per cell


### Output data
The code generates output .tsv files at different stages:
1. peptide buildup rates (log intensities, preptide ID split up in protein ID and peptide sequence, number of time points going into the fit, resulting rate)
2. tidied up peptide buildup rates (same as above, removed spurious rates smaller 1e-5)
3. protein buildup rates: peptide buildup rates averaged over (tidy) peptides corresponding to the same protein (swissprot id, rate, spread over which the different peptide estimates going into the rate range)
4. protein synthesis rates: buildup rates normalized by protein abundance per cell (swissprot id, HGNC gene symbo, ensembl gene id, protein copy number per cell averaged over replicates, protein synthesis rate, spread of the protein synthesis rate)

