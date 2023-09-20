import numpy as np 
import os
import sys
import pdb




frq_file = sys.argv[1]


f = open(frq_file)
head_count = 0

min_af = 1.0
for line in f:
	line = line.rstrip()
	data = line.split()
	if len(data) != 6:
		print('assumption erroror')
		pdb.set_trace()
	if head_count == 0:
		head_count = head_count + 1
		continue
	line_af = float(data[4])
	if line_af > .5:
		line_af = 1.0 - line_af
	if line_af < min_af:
		min_af = line_af

f.close()

print('minimum af is: ' + str(min_af))