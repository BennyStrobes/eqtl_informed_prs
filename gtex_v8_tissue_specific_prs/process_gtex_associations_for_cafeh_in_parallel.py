import numpy as np
import pdb
import sys
import os



def get_genes_on_chromosome(cafeh_gene_list, chrom_num):
    f = open(cafeh_gene_list)
    counter = 0
    head_count = 0
    genes = []
    gene_dicti = {}
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        counter = counter + 1
        if head_count == 0:
            head_count = head_count + 1
            continue
        gene_name = data[0]
        line_chrom_num = data[1]
        if line_chrom_num != 'chr' + chrom_num:
            continue
        genes.append(gene_name)
        gene_dicti[gene_name] = []
    f.close()
    return np.asarray(genes)


def get_list_of_genotyped_variants(variant_file):
    f = open(variant_file)
    dicti = {}
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        variant_id = data[0]
        dicti[variant_id] = 1
    f.close()
    return dicti

def array_to_dictionary_of_lists(ordered_genes):
    dicti = {}
    for gene in ordered_genes:
        dicti[gene] = []
    return dicti

def get_gtex_tissues(gtex_tissue_file):
    f = open(gtex_tissue_file)
    head_count = 0
    arr = []
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        arr.append(data[0])
    f.close()
    return np.asarray(arr)

def fill_in_gene_dicti_from_single_tissue_summary_stats(sum_stats_file, gene_dicti, gtex_tissue):
    f = open(sum_stats_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        gene_id = data[1]
        if gene_id not in gene_dicti:
            print('miss')
            pdb.set_trace()
            continue
        variant_id = data[0]
        beta = float(data[2])
        t_stat = float(data[3])
        beta_std_err = beta/t_stat
        pvalue = data[4]

        stringer = gtex_tissue + ',' + gene_id + ',' + variant_id + ',' + pvalue + ',' + str(beta) + ',' + str(beta_std_err)
        gene_dicti[gene_id].append(stringer)
    f.close()
    return gene_dicti

def fill_in_gene_dicti(gtex_tissues, chrom_num, eqtl_summary_stats_dir, gene_dicti):
    for gtex_tissue in gtex_tissues:
        print(gtex_tissue)
        sum_stats_file = eqtl_summary_stats_dir + gtex_tissue + '_chr' + chrom_num + '_matrix_eqtl_results.txt'

        gene_dicti = fill_in_gene_dicti_from_single_tissue_summary_stats(sum_stats_file, gene_dicti, gtex_tissue)
    return gene_dicti


cafeh_gene_list = sys.argv[1]
processed_gtex_associations_dir = sys.argv[2]
chrom_num = sys.argv[3]
eqtl_summary_stats_dir = sys.argv[4]
gtex_tissue_file = sys.argv[5]
genotype_reference_panel_dir = sys.argv[6]

# Extract ordered list of gtex tissues
gtex_tissues = get_gtex_tissues(gtex_tissue_file)

# Extract ordered list of genes on this chromosome
ordered_genes = get_genes_on_chromosome(cafeh_gene_list, chrom_num)


# Create dictionary where each key is a gene in ordered_genes and each value is an empty list
gene_dicti = array_to_dictionary_of_lists(ordered_genes)
# Fill in the dictionary with each elent in list is a string corresponding to info on a cis snp
gene_dicti = fill_in_gene_dicti(gtex_tissues, chrom_num, eqtl_summary_stats_dir, gene_dicti)


# Loop through genes
for gene_index, test_gene in enumerate(ordered_genes):
    # file containing variants for this gene
    variant_file = genotype_reference_panel_dir + 'genotype_reference_panel_' + test_gene + '.txt'
    genotyped_variant_list = get_list_of_genotyped_variants(variant_file)
    

    # print to output file
    output_file = processed_gtex_associations_dir + test_gene + '_associations.csv'
    t = open(output_file,'w')
    # Header
    t.write('tissue,gene_id,variant_id,pvalue,beta,beta_std_err\n')

    # Extract cis-snps for this gene
    gene_test_arr = gene_dicti[test_gene]

    # Print cis snps to output file
    for gene_test_string in gene_test_arr:
        test_variant_id = gene_test_string.split(',')[2]
        # Filter out variants with AF < THRESH across all gtex european samples
        # THRESH IS set to .005
        # Also filter out variants that are not fully observed across all gtex european samples (is about 1% of variants)
        if test_variant_id not in genotyped_variant_list:
            continue
        t.write(gene_test_string + '\n')
    t.close()

