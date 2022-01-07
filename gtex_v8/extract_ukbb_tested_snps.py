import numpy as np 
import os
import sys
import pdb
import gzip





sumstat_file = sys.argv[1]
output_root = sys.argv[2]

output_file_handles = {}
for chrom_num in range(1,23):
	chrom_string = str(chrom_num)
	output_file_handles[chrom_string] = open(output_root + chrom_string + '.txt','w')

f = gzip.open(sumstat_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	line_string = data[1]

	# Get allele frequency
	af = float(data[6])
	info = float(data[7])

	# Get maf from af
	if af > .5:
		maf = 1.0 - af
	else:
		maf = af
	if maf < .001 or info < .6:
		continue
	if line_string not in output_file_handles:
		print(line_string)
		continue
	snp_id = data[0]
	output_file_handles[line_string].write(snp_id + '\n')

f.close()




for chrom_num in range(1,23):
	chrom_string = str(chrom_num)
	output_file_handles[chrom_string].close()