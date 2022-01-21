import numpy as np 
import os
import sys
import pdb





genotype_file = sys.argv[1]
snp_loc_file = sys.argv[2]


f = open(genotype_file)
t = open(snp_loc_file,'w')
t.write('snp\tchr\tpos\n')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	variant_id = data[0]
	variant_info = variant_id.split('_')
	chrom_num = variant_info[0]
	pos = variant_info[1]
	t.write(variant_id + '\t' + chrom_num + '\t' + pos + '\n')
f.close()
t.close()
