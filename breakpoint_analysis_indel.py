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

INDEL_DIR = "/Users/vsubasri/research/indels_bed"
SNV_DIR = "/Users/vsubasri/research/SNVs_bed"
STATS = "/Users/vsubasri/research/snvs_breakpoint_10kb_indels_stats.tsv"
SUMMARY = "/Users/vsubasri/research/snvs_breakpoint_10kb_indels.tsv"
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
		df.to_csv(output, sep='\t', index=None)
	return df

def flanking_mutations(flank, indel_bed, indel_file, snv_dir):
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

def regex(search_term):
        regex = '(.*?)_annotated_filtered_indel.Bed'
        found = re.search(regex, search_term)
        if found:
                newname = found.group(1)  
                return newname

def distance_flanking(df):
	df['diff'] = abs(df[df.columns[5]] - df[df.columns[3]])
	return df

def compute_metrics(x):
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


#parser = argparse.ArgumentParser(description='Determine thresholds of breakpoint analysis')
#parser.add_argument('flank', metavar='N', type=int, nargs='+', help='distance flanking the breakpoint to measure snvs')
#parser.add_argument('--summary', dest='summary', action='store_const', const=summary:'summary', help='summary file', default='none')
#args = parser.parse_args()
