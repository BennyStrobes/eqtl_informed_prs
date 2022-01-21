import numpy as np 
import os
import sys
import pdb





genotype_stem = sys.argv[1]

vcf_file = genotype_stem + '.vcf'
dosage_file = genotype_stem + '_dosage.txt'



f = open(vcf_file)
t = open(dosage_file,'w')
for line in f:
	line = line.rstrip()
	data = line.split()
	if line.startswith('##'):
		continue
	if line.startswith('#CHROM'):
		# header
		raw_sample_names = data[9:]
		sample_names = []
		for raw_sample_name in raw_sample_names:
			sample_names.append(raw_sample_name.split('_')[1])
		sample_names = np.asarray(sample_names)
		t.write('id\t' + '\t'.join(sample_names) + '\n')
		continue
	variant_id = data[2]
	line_alt = data[4]
	variant_id_alt = variant_id.split('_')[3]
	if variant_id_alt == line_alt:
		flipped = False
	else:
		flipped = True

	genotypes = data[9:]
	if len(genotypes) != len(sample_names):
		print('assumption eroror')
		pdb.set_trace()
	dosages = []

	for genotype in genotypes:
		geno_a1 = genotype.split('/')[0]
		geno_a2 = genotype.split('/')[1]
		if geno_a1 != '0' and geno_a1 != '1':
			print('assumption eroror')
			pdb.set_trace()
		if geno_a2 != '0' and geno_a2 != '1':
			print('assumption erooro')
			pdb.set_trace()
		dosage = int(geno_a1) + int(geno_a2)
		dosages.append(dosage)
	dosages = np.asarray(dosages)
	if flipped == True:
		dosages = 2 - dosages
	t.write(variant_id + '\t' + '\t'.join(dosages.astype(str)) + '\n')
f.close()
t.close()