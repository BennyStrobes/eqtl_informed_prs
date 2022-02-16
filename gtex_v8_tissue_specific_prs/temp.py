import numpy as np 
import os
import sys
import pdb






gtex_tissue_file = sys.argv[1]
trait_name = sys.argv[2]
processed_bivariate_cafeh_input_dir = sys.argv[3]
bivariate_cafeh_output_dir = sys.argv[4]



gene_file = processed_bivariate_cafeh_input_dir + trait_name + '_processed_gene_list.txt'

t = open(bivariate_cafeh_output_dir + trait_name + '_num_components_comparison.txt','w')
t.write('full_cafeh\ttrait_only_cafeh\n')
f = open(gene_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene_name = data[0]
	file_name = bivariate_cafeh_output_dir + 'cafeh_debug_results_' + trait_name + '_' + gene_name + '_num_active_comparison.txt'
	if os.path.exists(file_name) == False:
		continue
	aa = np.loadtxt(file_name, dtype=str,delimiter='\t')
	t.write('\t'.join(aa[1,:]) + '\n')
f.close()
t.close()