#exec(open("GoodmanData.py").read())

from Bio import SeqIO
import csv
import pandas as pd
from pandas import *
from collections import defaultdict
import os
import sys
import glob
from pathlib import Path

data = read_csv("GoodmanDataset.csv")
CDS_seq = data['CDS.seq'].tolist()
Prot_FCC = data['Prot.FCC'].tolist()


codon_number_list=["CodonNumber"]
for i in range(len(CDS_seq)):
	for j in range(1,12):
		codon_number_list.append(j)

nuc3_list=["Nucleotide"]
for sequence in CDS_seq:
	nuc3 = sequence[2:33:3]
	for i in nuc3:
		nuc3_list.append(i)

Prot_FCC_list=["Prot.FCC"]
Prot_FCC_smalllist=[]
for i in Prot_FCC:
	for j in range(11):
		Prot_FCC_list.append(i)

with open('GoodmanDataTable.csv', 'w',newline='') as file:
	writer = csv.writer(file)
	for i in range(len(codon_number_list)):
		print(f"row{i} created")
		content = [codon_number_list[i],nuc3_list[i],Prot_FCC_list[i]]
		writer.writerow(content)
