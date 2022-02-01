import numpy as np 
import os
import sys
import pdb



def create_eqtl_test_info_dictionary(composit_tissues, num_composit_tissues, pseudotissue_eqtl_dir, chrom_num):
	dicti = {}
	for tissue in composit_tissues:
		tissue_chrom_eqtl_file = pseudotissue_eqtl_dir + tissue + '_chr' + chrom_num + '_matrix_eqtl_results.txt'
		f = open(tissue_chrom_eqtl_file)
		head_count = 0
		used_tests = {}
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			variant_id = data[0]
			gene_name = data[1]
			test_name = variant_id + ':' + gene_name
			if test_name in used_tests:
				print('assumption eroror')
				pdb.set_trace()
			used_tests[test_name] = 0
			beta = float(data[2])
			t_stat = float(data[3])
			beta_std_err = beta/t_stat
			if test_name not in dicti:
				dicti[test_name] = {}
				dicti[test_name]['beta'] = []
				dicti[test_name]['beta_std_err'] = []
			dicti[test_name]['beta'].append(beta)
			dicti[test_name]['beta_std_err'].append(beta_std_err)
		f.close()
	return dicti

def meta_analysis(effects, se, method='random', weights=None):
	# From Omer Weissbrod
	assert method in ['fixed', 'random']
	d = effects
	variances = se**2

	#compute random-effects variance tau2
	vwts = 1.0 / variances
	fixedsumm = vwts.dot(d) / vwts.sum()
	Q = np.sum(((d - fixedsumm)**2) / variances)
	df = len(d)-1
	tau2 = np.maximum(0, (Q-df) / (vwts.sum() - vwts.dot(vwts) / vwts.sum()))

	#defing weights
	if weights is None:
		if method == 'fixed':
			wt = 1.0 / variances
		else:
			wt = 1.0 / (variances + tau2)
	else:
		wt = weights

	#compute summtest
	summ = wt.dot(d) / wt.sum()
	if method == 'fixed':
		varsum = np.sum(wt*wt*variances) / (np.sum(wt)**2)
	else:
		varsum = np.sum(wt*wt*(variances+tau2)) / (np.sum(wt)**2)
	###summtest = summ / np.sqrt(varsum)

	summary=summ
	se_summary=np.sqrt(varsum)

	return summary, se_summary


pseudotissue_name = sys.argv[1]
composit_tissue_string = sys.argv[2]
pseudotissue_eqtl_dir = sys.argv[3]
chrom_num = sys.argv[4]


# Get names and number of composit tissues
composit_tissues = composit_tissue_string.split(',')
num_composit_tissues = len(composit_tissues)



# Create dictionary mapping test name to two arrays (one for effect sizes and one for standard errors)
eqtl_test_info = create_eqtl_test_info_dictionary(composit_tissues, num_composit_tissues, pseudotissue_eqtl_dir, chrom_num)


pseudotissue_chrom_meta_analysis_eqtl_file = pseudotissue_eqtl_dir + pseudotissue_name + '_meta_analyis_chr' + chrom_num + '_matrix_eqtl_results.txt'
t = open(pseudotissue_chrom_meta_analysis_eqtl_file,'w')
t.write('SNP\tgene\tbeta\tt-stat\n')

test_names = np.asarray([*eqtl_test_info])

for test_name in test_names:
	# Get test names
	test_info = test_name.split(':')
	variant_id = test_info[0]
	gene_id = test_info[1]
	# get vector of betas and standard errors corresponding to this test
	test_betas = np.asarray(eqtl_test_info[test_name]['beta'])
	test_beta_std_errs = np.asarray(eqtl_test_info[test_name]['beta_std_err'])
	# Error checking
	if len(test_betas) != len(test_beta_std_errs):
		print('assumption eororor')
		pdb.set_trace()
	# Throw out tests not fully observed across all composit tissues
	if len(test_betas) != num_composit_tissues:
		continue

	meta_beta, meta_beta_std_err = meta_analysis(test_betas, test_beta_std_errs, method='fixed')

	meta_t_stat = meta_beta/meta_beta_std_err

	t.write(variant_id + '\t' + gene_id + '\t' + str(meta_beta) + '\t' + str(meta_t_stat) + '\n')

t.close()




