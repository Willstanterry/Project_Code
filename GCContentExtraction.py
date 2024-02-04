
from collections import Counter
import xlsxwriter
import csv
import numpy as np
import pandas as pd
import numpy as np

filename = input("Please enter filename: ")
	
gene = []
coding = []
gene = ""
coding13 = []
coding13_1=[]
coding13_2=[]
coding13_3=[]

#Open file and find each independent gene and checks they are valid in length and characters
openfile = open(filename,'r')
code = openfile.readlines()[1:]
for line in code:
	if line.startswith('>'):
		gene_length = len(gene)
		if gene_length > 39 and gene_length % 3 == 0:
			coding.append(gene)
			gene = str("")

	else:
		for i in line:
			if i == "A" or i=="T" or i=="G" or i=="C":
				gene = gene + str(i)
coding.append(gene)



#Whole Genome
def codon_counter(gene):
	if len(gene) != 39:
		print(gene)
	return [len(gene[i:i+3].replace("A","").replace("T","")) for i in range(0, 39,3)]

for gene in coding:
    coding13.append(codon_counter(gene[:39]))
codons = np.array(coding13)
mean = codons.sum(axis=0)/len(codons)
df = pd.DataFrame(mean)
df.to_csv('GenomeData.csv')


#GC1
def codon_counter1(gene):
	return [len(gene[i].replace("A","").replace("T","")) for i in range(0, 39,3)]

for gene in coding:
	coding13_1.append(codon_counter1(gene[:39]))
codons1 = np.array(coding13_1)
mean1 = codons1.sum(axis=0)/len(codons1)
df = pd.DataFrame(mean1)
df.to_csv('GenomeDataGC1.csv')


#GC2
def codon_counter2(gene):
	return [len(gene[i+1].replace("A","").replace("T","")) for i in range(0, 39,3)]

for gene in coding:
	coding13_2.append(codon_counter2(gene[:39]))
codons2 = np.array(coding13_2)
mean2 = codons2.sum(axis=0)/len(codons2)
df = pd.DataFrame(mean2)
df.to_csv('GenomeDataGC2.csv')


#GC3
def codon_counter3(gene):
	return [len(gene[i+2].replace("A","").replace("T","")) for i in range(0, 39,3)]

for gene in coding:
	coding13_3.append(codon_counter3(gene[:39]))
codons3 = np.array(coding13_3)
mean3 = codons3.sum(axis=0)/len(codons3)
df = pd.DataFrame(mean3)
df.to_csv('GenomeDataGC3.csv')


#exec(open("GenomeExtraction.py").read())
