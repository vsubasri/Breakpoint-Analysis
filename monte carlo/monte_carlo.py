import argparse 
import random
import pandas as pd
import mutation_rate
import numpy as np
import os


##things to consider: read depth, GC count 

DIR = "/Users/vsubasri/Medical Biophysics 2017:2018/Rotation Project 3/SNVs_bed"
OUTPUT = "/Users/vsubasri/Medical Biophysics 2017:2018/Rotation Project 3/monte_carlo/"

N = 1501753	#3137161264
S = 10

def decode_base(base):
	decoded_base = ''
	if base == 0: 
		decoded_base = 'A' 
	elif base == 1: 
		decoded_base = 'G'
	elif base == 2: 
		decoded_base = 'C'
	else: 
		decoded_base = 'T'
	return decoded_base 


def decode_bases(sequence):
	decoded_sequence = []
	for base in range(len(sequence)): 
		base = decode_base(sequence[base])
		decoded_sequence.append(base)
	return decoded_sequence

##Generate random sequence of bases for flanking regions 
def generate_sequence(): 
	sequence = []
	for i in range(N): 
		base = random.randint(0,4)
		sequence.append(base)
	decoded_sequence = decode_bases(sequence)
	return decoded_sequence 

##Run simulation of mutations in N bp 

def random_simulation(file):
	sequence = generate_sequence()
	mutated = []
	snvs = []
	rate = mutation_rate.calculate_rate(file)
	for i in range(len(sequence)):
		snv =  np.random.choice([0, 1], 1, p=[(1-rate), rate])			##1 represents mutation, 0 represents no mutation 
		snvs.append(snv)
		if snv == 1: 
			base = random.randint(0,4)
			mutated.append(decode_base(base))
		if snv == 0:
			mutated.append(sequence[i])
	count = [sequence, mutated]
	return count

def run_simulations(file):
	simulations = []
	for i in range(S):
		simulation = random_simulation(file)
		simulations.append(simulation)
	return simulations

def combine_mutations(file):
	combined = []
	simulations = run_simulations(file)
	for s in range(len(simulations)):
		pre = simulations[s][0]
		post = simulations[s][1]
		entry = []
		for i in range(len(pre)): 
			entry.append(pre[i] + post[i])
		combined.append(entry)
	return combined

def mutations(file):
	combined = combine_mutations(file)
	df = pd.DataFrame(index=list(range(0, len(combined))), columns=['AT', 'AG', 'AC', 'AA', 'TA', 'TC', 'TG', 'TT', 'CT', 'CA', 'CG', 'CC', 'GT', 'GA', 'GC', 'GG'])
	df.fillna(0, inplace=True)
	for i in range(len(combined)):
		for base in combined[i]:
			df.loc[i][base] += 1
	df.to_csv(OUTPUT+mutation_rate.regex(file)+'.tsv', sep='\t', index=None)
	return df

def compute_statistics(df):
	count = df[['AT', 'AG', 'AC', 'AA', 'TA', 'TC', 'TG', 'TT', 'CT', 'CA', 'CG', 'CC', 'GT', 'GA', 'GC', 'GG']].count()
	print(count)

def simulations_all():
	simulations = []
	files = next(os.walk(DIR))[2]
	for i in range(len(files)):
		df = mutations(files[i])
		compute_statistics(df)
		break

simulations_all()
