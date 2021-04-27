#!/usr/bin/env python3

#Author: Devang Haresh Liya
#Cite: 
#Website: https://qpromoters.com

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import motifs
import matplotlib.pyplot as plt
import numpy as np

print("   ____  _____                           _                ")
print("  / __ \\|  __ \\                         | |               ")
print(" | |  | | |__) | __ ___  _ __ ___   ___ | |_ ___ _ __ ___ ")
print(" | |  | |  ___/ '__/ _ \\| '_ ` _ \\ / _ \\| __/ _ \\ '__/ __|")
print(" | |__| | |   | | | (_) | | | | | | (_) | ||  __/ |  \\__ \\")
print("  \\___\\_\\_|   |_|  \\___/|_| |_| |_|\\___/ \\__\\___|_|  |___/")
print("                                                     ")
                                                          
print("Your go-to tool for predicting the strength of Saccharomyces cerevisiae promoters!\n                                                          ")

def calculate_score(sequence, pssm):
	score = 0
	for i, nt in enumerate(sequence):
		score = score + pssm[nt, i]
	return score

path = 'EPDPromoters/' 

fn_all = 'all_promters_49_10.fa'	#All promoters from S_cerevisiae_epdnew database
file_in_all = path + fn_all
all_seq = {}	#list of all sequences from EPD database

with open(file_in_all, 'r') as f:
	for seq_record in SeqIO.parse(f, 'fasta'):
		pname = ''.join(seq_record.description.split()[1])
		all_seq[pname] = Seq(str(seq_record.seq).upper())

#Seperate dictionary for characterized promoters
characterized_names = ['TDH3_1', 'CCW12_1', 'PGK1_1', 'HHF2_1', 'TEF1_1', 'TEF2_1', 'HHF1_1', 'HTB2_1', 'RPL18B_1', 'ALD6_1', 'PAB1_1',
'RET2_1', 'RNR1_1', 'SAC6_1', 'RNR2_1', 'POP6_1', 'RAD27_1', 'PSP2_1']
characterized_seq = {pname: all_seq[pname] for pname in characterized_names}

#Create motif
moti = motifs.create(list(all_seq.values()))
pwm = moti.counts.normalize()
pssm = pwm.log_odds()

#Calculate scores
TDH3score = calculate_score(characterized_seq['TDH3_1'], pssm) #Strongest promoter in Lee et al. Will normalize by dividing with this score

in_epd = input("Is your promoter in Eukaryotic Promoter Databse (EPD)?[y/n] ").lower()
if in_epd == 'y':
	userpromoter = input("Please enter the EPDnew ID of your promoter: ").upper()
	while userpromoter not in all_seq.keys():
		userpromoter = input("Oops! Didn't find that in EPD.\nPlease enter the EPDnew ID of your promoter: ").upper()
	userseq = all_seq[userpromoter]
else:
	userseq = input("Please enter a 60 nucleotide sequence: ").upper()
	while (len(userseq) != 60) or (not all(i in 'ATCG' for i in userseq)):
		userseq = input("Oops! Looks like you entered an invalid sequence.\nMake sure your sequence is exactly 60 nucleotides long and it contains only letter A, T, C, G.\nPlease enter a 60 nucleotide sequence: ").upper()

userscore = calculate_score(userseq, pssm)
print("\n-----------------------Results-----------------------")
print("Your promoter sequence was: " + userseq)
print("Raw promoter score:", round(userscore,3))
print("Promoter score normalized by TDH3:", round(userscore/TDH3score, 3))
print("Predicted promoter strength (Liya et al):", round(0.93*(userscore/TDH3score) - 0.07, 3))


#Calculate scores of all EPD promoters
all_scores = {}
for pname, pseq in all_seq.items():
	all_scores[pname] = calculate_score(pseq, pssm)

#Calculate the scores of all characterized promoters
sorted_names = ['RET2_1', 'PAB1_1', 'POP6_1', 'RNR2_1', 'RAD27_1', 'PSP2_1', 'RNR1_1', 'SAC6_1', 'TEF2_1', 'HTB2_1', 'HHF1_1', 'PGK1_1',
'TEF1_1', 'HHF2_1', 'CCW12_1', 'RPL18B_1', 'ALD6_1', 'TDH3_1']
characterized_scores = []
for pname in sorted_names:
	characterized_scores.append(all_scores[pname])

plt.figure(figsize=(10,4))

plt.subplot(1,2,1)
plt.scatter(range(len(characterized_scores)), characterized_scores, color='blue', label=r"Lee $et\:al.$ promoters")
plt.xticks(range(len(characterized_scores)), sorted_names, fontsize='x-small', rotation=45)

plt.axhline(userscore, color='red', label="Your promoter")
plt.ylabel("-49 to 10 score")
plt.title(r"Comparison with Lee $et\:al.$ promoters")
plt.legend()

plt.subplot(1,2,2)
# find out in which bin belongs the position where you want the label
ybins, xbins, _ = plt.hist(all_scores.values(), bins=50)
labeled_bin = userscore
ind_bin = np.where(xbins >= labeled_bin)[0]
if len(ind_bin) > 0 and ind_bin[0] > 0:
	# get position and value of the bin
	x_bin = xbins[ind_bin[0]-1]/2. + xbins[ind_bin[0]]/2.
	y_bin = ybins[ind_bin[0]-1]
	# add the arrow
	plt.annotate("Your promoter", xy=(x_bin, y_bin + 20), xycoords='data', xytext=(x_bin, y_bin + 200), textcoords='data',
		arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))

plt.xlabel("Promoter scores")
plt.title(r"Comparison with 5117 $S. cerevisiae$ promoters")
plt.tight_layout()
plt.show()

print("-----------------------------------------------------")