import argparse
import csv
import pandas as pd
import os
import subprocess
import re 
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

sns.set(color_codes=True)

INDEL_DIR = "/Users/vsubasri/research/filtered_tab_files_bed"
SNV_DIR = "/Users/vsubasri/research/SNVs_bed"
OUTPUTDIR="/Users/vsubasri/research/breakpoint_analysis/output/"
STATS = os.path.join(OUTPUTDIR,"snvs_breakpoint_10kb_stats.tsv")
SUMMARY = os.path.join(OUTPUTDIR,"snvs_breakpoint_10kb.tsv")
DISTRIBUTION= os.path.join(OUTPUTDIR,"snvs_breakpoint_10kb_distribution.tsv")
HIGH=os.path.join(OUTPUTDIR,"snvs_breakpoint_10kb_high.tsv")

FLANK = 10000

def process_mutations(indel_dir, flank, output):
	files = next(os.walk(indel_dir))[2]
	snvs = []
	for file in files:
		bed = pd.read_table(indel_dir+"/"+file, sep='\t', index_col=None)
		mutations = flanking_mutations(flank, bed, file, SNV_DIR)
		snvs.append(mutations)
	df = pd.concat(snvs)
	df.columns = ['sample', 'chr', 'start', 'end', 'score', 'breakpoint', 'diff']
	if output == "none":
		pass
	else:
##		df.to_csv(output, sep='\t', index=None)
	return df

def flanking_mutations_indel(flank, indel_bed, indel_file, snv_dir):
	snvs = []
	print(indel_file.split("_indel")[0])
	for indel in indel_bed.index:
		chr = indel_bed['annovar_chr'][indel]	##indel breakpoint chromosome
		indel_breakpoint = indel_bed['annovar_start'][indel]	##indel breakpoint location
		snv_file = indel_file.split("_indel")[0] + "_clipped" + indel_file.split("_indel")[1]
		snv_bed = pd.read_table(snv_dir+"/"+snv_file, sep='\t', index_col=None)
		for snv in snv_bed.values:
			score = snv[3] + snv[4]
			snv_entry = [regex(indel_file), snv[0], snv[1], snv[2], score, indel_breakpoint]
			if (snv_entry[1] == chr):
				plus = snv_entry[2] + flank
				minus = snv_entry[2] - flank
				if (plus > indel_breakpoint and minus < indel_breakpoint):
					snvs.append(snv_entry)#	print(snvs)
	df = distance_flanking(pd.DataFrame(snvs))
	return df

def flanking_mutations(flank, indel_bed, indel_file, snv_dir):
	snvs = []
	for indel in indel_bed.index:
		chr1 = indel_bed['Chromosome1'][indel]	##indel breakpoint chromosome
		chr2 = indel_bed['Chromosome2'][indel]
		indel_breakpoint1 = indel_bed['Position1'][indel]	##indel breakpoint location 1
		indel_breakpoint2 = indel_bed['Position2'][indel] ##indel breakpoint location 2
		snv_file = indel_file.split(".")[0] + "_annotated_filtered_clipped.Bed" ##find snvs for corresponding sample
		snv_bed = pd.read_table(snv_dir+"/"+snv_file, sep='\t', index_col=None)
		for snv in snv_bed.values:
			score = snv[3] + snv[4] ##type of mutation ie. A > G, C > T, etc.
			snv_entry1 = [regex(indel_file), snv[0], snv[1], snv[2], score, indel_breakpoint1]
			snv_entry2 = [regex(indel_file), snv[0], snv[1], snv[2], score, indel_breakpoint2]
			##if snv chromosome is the same as breakpoint1 chromosome
			if (snv_entry1[1] == regex(chr1)):
				plus_bp1 = snv_entry1[2] + flank
				minus_bp1 = snv_entry1[2] - flank
				if (plus_bp1 > indel_breakpoint1 and minus_bp1 < indel_breakpoint1):		
					snvs.append(snv_entry1)#	print(snvs)
			##if snv chromosome is the same as breakpoint2 chromosome		
			elif (snv_entry2[1] == regex(chr2)):
				plus_bp2 = snv_entry2[2] + flank
				minus_bp2 = snv_entry2[2] - flank
				if (plus_bp2 > indel_breakpoint2 and minus_bp2 < indel_breakpoint2):	
					snvs.append(snv_entry2)#	print(snvs)
	if not snvs:
		df=None
		print("No SNVs flanking " + os.path.basename(indel_file))
	else:
		print("SNVs found flanking breakpoints in " + indel_file.split("_indel")[0])
		df = distance_flanking(pd.DataFrame(snvs))
		return df

def regex(search_term):
	if ("_annotated_filtered_indel.Bed" in search_term):
		regex = '(.*?)_annotated_filtered_indel.Bed'
	elif ("filtered.Bed" in search_term):
		regex = '(.*?).(inv|tra|del|dup).filtered.Bed'
	else:
		regex = 'chr(.*)'
	found = re.search(regex, search_term)
	if found:
		newname = found.group(1) 
		return newname

def distance_flanking(df):
	df['diff'] = abs(df[df.columns[5]] - df[df.columns[3]])
	return df

def compute_metrics(x):
	distribution = pd.DataFrame(x.groupby(['sample','breakpoint'])['score'].count())
	distribution.to_csv(DISTRIBUTION, sep='\t')
	high = distribution[distribution['score']>2]
	high.to_csv(HIGH, sep='\t')
	print(high)
	####
	metrics = x.groupby('sample')['diff'].agg(['count','mean', 'std'])
	close = x[x['diff'] < 5].groupby('sample')['diff'].count()
	far = x[x['diff'] > 5].groupby('sample')['diff'].count()
	TA = x[x['score'].str.contains('TA|AT')].groupby('sample')['score'].count()
	TG = x[x['score'].str.contains('TG|AC')].groupby('sample')['score'].count()
	TC = x[x['score'].str.contains('TC|AG')].groupby('sample')['score'].count()
	CT = x[x['score'].str.contains('CT|GA')].groupby('sample')['score'].count()
	CA = x[x['score'].str.contains('CA|GT')].groupby('sample')['score'].count()
	CG = x[x['score'].str.contains('CG|GC')].groupby('sample')['score'].count()
	df = pd.DataFrame(pd.concat([metrics, close, far, TA, TG, TC, CT, CA, CG], axis=1))
	df.columns=['count', 'mean', 'std', 'close(<5bp)', 'far(>5bp)', 'TA', 'TG', 'TC', 'CT', 'CA', 'CG']
	return df

def summary_stats(indel_dir, flank, output_stats, output_summ):
#	df = pd.read_csv(SUMMARY, sep='\t')
	df = process_mutations(indel_dir, flank, output_summ) 
	concat = compute_metrics(df)
	concat.to_csv(output_stats, sep='\t')
	plot_count(concat.iloc[0:, 3:])

def plot_count(concat):
	concat['sample'] = concat.index.tolist()
	count = pd.melt(concat.iloc[0:,2:], id_vars="sample", var_name="measure", value_name="count")
	sns_plot = sns.factorplot(x='sample', y='count', hue='measure', data=count, kind='bar', size=15, orient="v", legend=False)
	sns_plot.set_xticklabels(rotation=30)
	plt.legend(bbox_to_anchor=(1.01, 0.9), loc=2, borderaxespad=0.5)
	sns_plot.savefig("output.png")


summary_stats(INDEL_DIR, FLANK, STATS, SUMMARY)	

