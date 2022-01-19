



ukbb_sample_names_dir="$1"
ukbb_download_data="$2"
ukbb_sumstats_dir="$3"
processed_ukbb_genotype_dir="$4"


ukbb_test_snps_output_root=$processed_ukbb_genotype_dir"ukbb_tested_snps_chr_"
sumstat_file=$ukbb_sumstats_dir"bolt_337K_unrelStringentBrit_MAF0.001_v3.blood_WHITE_COUNT.bgen.stats.gz"
if false; then
python extract_ukbb_tested_snps.py $sumstat_file $ukbb_test_snps_output_root
fi


sample_names_keep_file=${ukbb_sample_names_dir}keep.non_british_european_122K.txt
chrom="22"
snp_names_keep_file=${ukbb_test_snps_output_root}${chrom}".txt"

plink2 --bgen ${ukbb_download_data}ukb_imp_chr${chrom}_v3.bgen 'ref-last' --sample ${ukbb_download_data}ukb1404_imp_chr1_v3_withdrawn3.sample --keep ${sample_names_keep_file} --extract ${snp_names_keep_file} --threads 1 --make-pgen --out ${processed_ukbb_genotype_dir}plink2_chr${chrom}.UKB500K