from Bio import SeqIO
import csv
import pandas as pd
from pandas import *
from collections import defaultdict
import os
import sys
import glob
from pathlib import Path

path = os.getcwd()
srch_path = os.path.join(path,"orthologs_foreve", "output_sequences")
files_aln = glob.glob(os.path.join(srch_path,"*.fta_aln"))

changes_ab = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
file_counter = 0

for file in files_aln:
	file_counter +=1
	open(file,'r')
	sequences = []
	
	for seq_record in SeqIO.parse(file, 'fasta'):
		sequences.append(seq_record.seq)
	
	if len(sequences)!=2:
		file_counter = file_counter - 1
		continue
	
	if (len(sequences[0]) + len(sequences[1])) > 120:
		for j in range(0,20):
			i = (j*3)+2
			if sequences[0][i] != sequences[1][i]:
				changes_ab[j]+=1
	else:
		file_counter = file_counter - 1
		
conservation_ab_perc_list = ["Percentage"]

for change in changes_ab:
	change_ab_perc = (change/file_counter)*100
	conservation_ab_perc = 100 - change_ab_perc
	print(f"There is a {conservation_ab_perc}% absolute convservation")
	conservation_ab_perc_list.append(conservation_ab_perc)

Codon_Number = ["CodonNumber",1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

with open('Absolute_Conservation_Percentage.csv', 'w', newline='') as csvfile:
	writer = csv.writer(csvfile)
	for i in range((len(Codon_Number)-1)):
		writer.writerow([Codon_Number[i],conservation_ab_perc_list[i]])

#############

path = os.getcwd()
srch_path = os.path.join(path,"orthologs_foreve", "output_sequences")
files_aln = glob.glob(os.path.join(srch_path,"*.fta_aln"))

changes_ATpres = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
changes_GCpres = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
file_counter = 0

for file in files_aln:
	file_counter +=1
	open(file,'r')
	sequences = []
	
	for seq_record in SeqIO.parse(file, 'fasta'):
		sequences.append(seq_record.seq)
	
	if len(sequences)!=2:
		file_counter = file_counter - 1
		continue
	
	sq1=[]
	sq2=[]	
	
	if (len(sequences[0]) + len(sequences[1])) > 120:
		for j in range(0,20):
			i = (j*3)+2
			if sequences[0][i]!=sequences[1][i]:
			
				if sequences[0][i] == "a":
					if sequences[1][i] != "t":
						changes_ATpres[j]+=1
				elif sequences[0][i] == "t":
					if sequences[1][i] != "a":
						changes_ATpres[j]+=1
				elif sequences[0][i] == "c":
					if sequences[1][i] != "g":
						changes_GCpres[j]+=1
				elif sequences[0][i] == "g":
					if sequences[1][i] != "c":
						changes_GCpres[j]+=1
				else:
					continue
	else:
		file_counter = file_counter - 1


conservation_ATpres_perc_list = ["Percentage"]
for change in changes_ATpres:
	change_ATpres_perc = (change/file_counter)*100
	conservation_ATpres_perc = 100 - change_ATpres_perc
	print(f"There is a {conservation_ATpres_perc}% AT-group preservation conservation")
	conservation_ATpres_perc_list.append(conservation_ATpres_perc)

Codon_Number = ["CodonNumber",1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

with open('Preserving_AT_Conservation_Percentage.csv', 'w', newline='') as csvfile:
	writer = csv.writer(csvfile)
	for i in range((len(Codon_Number)-1)):
		writer.writerow([Codon_Number[i],conservation_ATpres_perc_list[i]])


conservation_GCpres_perc_list = ["Percentage"]
for change in changes_GCpres:
	change_GCpres_perc = (change/file_counter)*100
	conservation_GCpres_perc = 100 - change_GCpres_perc
	print(f"There is a {conservation_GCpres_perc}% GC-group preservation conservation")
	conservation_GCpres_perc_list.append(conservation_GCpres_perc)
	
Codon_Number = ["CodonNumber",1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

with open('Preserving_GC_Conservation_Percentage.csv', 'w', newline='') as csvfile:
	writer = csv.writer(csvfile)
	for i in range((len(Codon_Number)-1)):
		writer.writerow([Codon_Number[i],conservation_GCpres_perc_list[i]])

