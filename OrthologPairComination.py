#exec(open("OrthologData.py").read())


from Bio import SeqIO
import csv
import pandas as pd
from pandas import *
from collections import defaultdict
import os
import sys

#Creates lists of identifying words for gene
data = read_csv("orthologs_foreve/3528-Proteins-from-a-list-2.csv")
Gene_Accession1 = data['Gene Accession1'].tolist()
Gene_Accession2 = data['Gene Accession2'].tolist()


#Other One Search
#Has SDY_RS123 number supplied and appends description and sequence to list for EC only
filenameIn2 = "/Volumes/ProjectCode/projectcode/Task2_OrthologPairs/orthologs_foreve/sequenceSD.txt"
sequences = [i for i in SeqIO.parse(filenameIn2,'fasta')]

sequence_desc = ""
wanted_sequence = []
filenumber = 0

for i in range(len(Gene_Accession1)):
	key2 = str(Gene_Accession1[i])
	for sequence in sequences:	
		sequence_desc = sequence.description
		if key2 in sequence_desc:				
			if key2 == "nan":
				filenumber = filenumber + 1
				continue
			found_sequence = sequence
			filename = f"/Volumes/ProjectCode/projectcode/Task2_OrthologPairs/orthologs_foreve/output_sequences/file{i}.fta"
			with open(filename, 'w') as file:
				file.write(">" + str(found_sequence.id)+ "\n")
				file.write(str(found_sequence.seq))
				print(f"Created file {i}")
				file.close()
		else:
			sequence_desc = ""



#E. Coli search
#Has b0000 number supplied and appends description and sequence to list for EC only
filenameIn1 = "/Volumes/ProjectCode/projectcode/Task2_OrthologPairs/orthologs_foreve/sequenceEC.txt"
sequences = [i for i in SeqIO.parse(filenameIn1,'fasta')]

sequence_desc = ""
wanted_sequence = []
filenumber = 0

for i in range(len(Gene_Accession2)):
	key1 = str(Gene_Accession2[i])
	for sequence in sequences:
		sequence_desc = sequence.description
		if key1 in sequence_desc:				
			if key1 == "nan":
				filenumber = filenumber + 1
				continue
			found_sequence = sequence
			filename = f"/Volumes/ProjectCode/projectcode/Task2_OrthologPairs/orthologs_foreve/output_sequences/file{i}.fta"
			with open(filename, 'a') as file:
				file.write("\n")
				file.write(">" + str(found_sequence.id)+ "\n")
				file.write(str(found_sequence.seq))
				print(f"Edited file {i}")
				file.close()
		else:	
			sequence_desc = ""	
		
path = os.getcwd()
srch_path = os.path.join(path,"orthologs_foreve", "output_sequences")
files_aln = glob.glob(os.path.join(srch_path,'*.fta'))

for fl in files_aln:
	print(f"Aligning {fl}")
	output_mafft = f"{fl}_aln"
	os.system('mafft --quiet %s > %s'%(fl,output_mafft))
	os.remove(fl)
	
	
	
	
