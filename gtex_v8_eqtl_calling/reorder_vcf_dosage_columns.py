import numpy as np 
import os
import sys
import pdb





file_stem = sys.argv[1]
indi_names_file = sys.argv[2]
sample_names_file = sys.argv[3]


sample_names = np.loadtxt(sample_names_file,dtype=str,delimiter='\t')

input_file = file_stem + '_dosage.txt'
output_file = file_stem + '_dosage_sample_ordered.txt'


head_count = 0

f = open(input_file)
t = open(output_file,'w')
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		file_names = np.asarray(data[1:])
		dicti = {}
		for index, file_name in enumerate(file_names):
			dicti[file_name] = index
		mapping_arr = []
		ordered_indis = []
		for sample_name in sample_names:
			mapping_arr.append(dicti[sample_name.split(':')[1]])
			ordered_indis.append(sample_name.split(':')[1])
		mapping_arr = np.asarray(mapping_arr)
		ordered_indis = np.asarray(ordered_indis)

		if np.array_equal(ordered_indis, file_names[mapping_arr]) == False:
			print('assumption eroror')
			pdb.set_trace()

		t.write(data[0] + '\t' + '\t'.join(sample_names) + '\n')
		continue
	variant_id = data[0]
	dosages = np.asarray(data[1:])
	reordered_dosages = dosages[mapping_arr]
	t.write(variant_id + '\t' + '\t'.join(reordered_dosages) + '\n')
f.close()
t.close()