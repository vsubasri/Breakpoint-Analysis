import pandas as pd
import os 
import subprocess
import re
import csv
from collections import OrderedDict

DIR = "/Users/vsubasri/Medical Biophysics 2017:2018/Rotation Project 3/SNVs_bed"
GENOME = float(3137161264)

OUTPUT = "/Users/vsubasri/Medical Biophysics 2017:2018/Rotation Project 3/mutation_rates.tsv"

def count_mutations(fname):
	os.chdir(DIR)
	cmd = ["wc", "-l", fname]
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	j = proc.communicate()[0]
	j = j.decode('utf-8').strip().split(' ')[0]
	return int(j)

def regex(search_term):
        regex = '(.*?)_annotated_filtered_clipped.Bed'
        found = re.search(regex, search_term)
        if found:
                newname = found.group(1)  
                return newname

def count_snvs(bed_dir):
	snvs = {}
	files = next(os.walk(bed_dir))[2]
	for file in files:
		snvs.update({file: float(count_mutations(file))/GENOME})
	return snvs

def calculate_rate(file):
	rate = float(count_mutations(file))/GENOME
	return rate

def output_rates(output_file):
	snvs = OrderedDict()
	snvs.update({'sample': 'mutation rate'})
	files = next(os.walk(DIR))[2]
	for file in files:
		snvs.update({regex(file): float(count_mutations(file))/float(genome)})
	with open(OUTPUT, 'wb') as f:
		writer = csv.writer(f, delimiter ='\t')
		for row in snvs.iteritems():
			writer.writerow(row)

output_rates(OUTPUT)
