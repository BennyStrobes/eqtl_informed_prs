import numpy as np 
import os
import sys
import pdb



input_bim = sys.argv[1]
output_bim = sys.argv[2]



f = open(input_bim)
t = open(output_bim, 'w')

for line in f:
	line = line.rstrip()
	data = line.split('\t')
	chrom_num = data[0]
	chrom_pos = data[3]
	a1 = data[4]
	a2 = data[5]
	gtex_id1 = 'chr' + chrom_num + '_' + chrom_pos + '_' + a1 + '_' + a2 + '_b38'
	gtex_id2 = 'chr' + chrom_num + '_' + chrom_pos + '_' + a2 + '_' + a1 + '_b38'
	#t.write(data[0] + '\t' + gtex_id1 + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5] + '\n')
	#t.write(data[0] + '\t' + gtex_id2 + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5] + '\n')
	t.write(gtex_id1 + '\n')
	t.write(gtex_id2 + '\n')
f.close()
t.close()