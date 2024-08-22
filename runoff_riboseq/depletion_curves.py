#!/usr/bin/env python

import numpy as np
import pandas as pd
import pysam
import json
import matplotlib.pyplot as plt



#file paths
tsv_file=("/scicore/home/zavolan/schlus0000/riboseq_pipeline/snakemake/prepare_annotation/results/homo_sapiens/transcript_id_gene_id_CDS.tsv")

time_files_dict = { 0 : ("trplct_10_Huh7_CHXHRT/transcripts.mapped.unique.sorted.sam","trplct_10_Huh7_CHXHRT/p_site_offsets/alignment_offset.json"), \
	80 : ("trplct_10_Huh7_HRT_1_3min/transcripts.mapped.unique.sorted.sam","trplct_10_Huh7_HRT_1_3min/p_site_offsets/alignment_offset.json"), \
	90 : ("trplct_10_Huh7_HRT_1_5min/transcripts.mapped.unique.sorted.sam","trplct_10_Huh7_HRT_1_5min/p_site_offsets/alignment_offset.json"), \
	100 : ("trplct_10_Huh7_HRT_1_6min/transcripts.mapped.unique.sorted.sam","trplct_10_Huh7_HRT_1_6min/p_site_offsets/alignment_offset.json"), \
	120 : ("trplct_10_Huh7_HRT_2min/transcripts.mapped.unique.sorted.sam","trplct_10_Huh7_HRT_2min/p_site_offsets/alignment_offset.json"), \
	150 : ("trplct_10_Huh7_HRT_2_5min/transcripts.mapped.unique.sorted.sam","trplct_10_Huh7_HRT_2_5min/p_site_offsets/alignment_offset.json"), \
	180 : ("trplct_10_Huh7_HRT_3min/transcripts.mapped.unique.sorted.sam","trplct_10_Huh7_HRT_3min/p_site_offsets/alignment_offset.json") }

output_filename = "output.txt"

#parameters
minlen = 3000
display_length = minlen
SL_binsize = 15
normwindow = 600

#define global dictionaries
CDS_start_ind = {}
CDS_stop_ind = {}
CDS_lengths = {}



'''
method reads the .tsv file and initializes: 
 - a dictionary with start locations of the CDS
 - a dictionary with end locations of the CDS
 - a dictionary with CDS lengths
'''
def read_tsv(tsv_file):
	with open(tsv_file) as CDS_coordinates:
		for line in CDS_coordinates:
			sp_line = line.strip().split("\t")
			transcript_ID = sp_line[0]
			CDS_start_ind[transcript_ID] = int(sp_line[2]) - 1
			CDS_stop_ind[transcript_ID] = int(sp_line[3])
			CDS_lengths[transcript_ID] = int(sp_line[3]) - int(sp_line[2]) + 1



'''
method reads sam/bam file
gives binned output curve
'''
def bam2meta_curve(sam_filename,json_filename,display_length=3000,min_CDS_length=3000,binsize=1):
	sam = pysam.AlignmentFile(sam_filename, 'rb')
	with open(json_filename, 'r') as json_file:
		 json_data = json_file.read()
	offset_obj = json.loads(json_data)

	counts = np.zeros(display_length//binsize)
	for read in sam.fetch():
		if read.is_reverse:
			continue
		
		transcript = read.reference_name
		read_length = len(read.seq)
		
		if CDS_lengths[transcript] >= min_CDS_length and str(read_length) in offset_obj:
			p_site_pos = (read.reference_start + offset_obj[str(read_length)] - CDS_start_ind[transcript])

			if p_site_pos >= 0 and p_site_pos < min(CDS_lengths[transcript],display_length):
				counts[p_site_pos//binsize] += 1
	
	return counts



'''
BEGIN SCRIPT
'''
read_tsv(tsv_file)


fig,ax = plt.subplots()

x_data = np.arange(0,display_length//SL_binsize)

w = open(output_filename,"w")

#print x array
for x in x_data:
	w.write("%d\t"%(x))
w.write("\n")

for t,(sfile,jsonfile) in time_files_dict.items():
	y_data = bam2meta_curve(sfile,jsonfile,display_length=display_length,min_CDS_length=minlen,binsize=SL_binsize)

	norm = sum(y_data[(minlen - normwindow)//SL_binsize:(minlen//SL_binsize)]) / (normwindow//SL_binsize)

	y_data = y_data / norm

	ax.plot(x_data,y_data,label=str(t)+"s")

	for y in y_data:
		w.write("%d\t" % (y))
	w.write("\n")

w.close()

plt.xticks([0,33.333,66.667,100,133.333,166.667,200],[0,500,1000,1500,2000,2500,3000])

plt.xlim([0,display_length//SL_binsize])
plt.ylim([0,3])

plt.xlabel('p-site position [nts]')
plt.ylabel('fraction of total counts')

plt.legend(loc="upper right")

plt.savefig('depletion_plot.pdf', bbox_inches='tight')

