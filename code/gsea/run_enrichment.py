from collections import defaultdict,OrderedDict
import numpy as np
import sys
from scipy.stats import fisher_exact
import os
from gsea import readgenesets, gsea
file_name = sys.argv[1] # File with test genes (one gene per line)
background_file_name = sys.argv[2] # File with background genes (one gene per line)
geneset_name = sys.argv[3] # filepath+filename of a geneset downloaded from msigdb. Need not parse msigdb file. Script takes care of parsing
save_name = sys.argv[4] # filepath+name for output file

""" Read test genes"""
fh = open(file_name,'r')
lines = fh.readlines()
fh.close()
test_genes = []
for line in lines:
	test_genes.append(line.strip('\n'))


""" Read background genes"""
fh = open(background_file_name,'r')
lines = fh.readlines()
fh.close()
all_genes = []
for line in lines:
	all_genes.append(line.strip('\n'))

		
all_genes = set(all_genes)
test_genes = set(test_genes)



test_out = gsea(test_genes, all_genes, geneset_name)

""" Save result to save_name"""
f = open(save_name,'w')
print("TestFile:", file_name, "\nBackgroundFile:", background_file_name, "\nGeneset:", geneset_name, "SaveFile:", save_name, file=f)
print("\ttest_inset\ttest_notinset\tbackground_inset\tbackground_notinset\toddsratio\tpvalue\tbonferroni-adjusted\ttestgenes_in_set", file=f)
for i in test_out:
	print(i,'\t','\t'.join(str(j) for j in test_out[i]), file=f)
f.close()
