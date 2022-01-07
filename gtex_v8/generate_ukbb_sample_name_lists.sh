


ukbb_download_data="$1"
ukbb_sample_names_dir="$2"



# make European-non-British sample list:
echo -e "\t Making non-British European remove list\n"
cat ${ukbb_download_data}/bolt.in_plink_but_not_imputed.FID_IID.976.txt \
    ${ukbb_download_data}/../sampleQC/remove.nonWhite.FID_IID.txt \
    ${ukbb_download_data}/w14048_20181016.FID_IID.txt \
    ${ukbb_download_data}/../sampleQC/samples_337K.txt |sort > ${ukbb_sample_names_dir}/remove.non_british_european.txt

echo -e "\t Sorting all samples\n"
cut -f 1,2 -d " " ${ukbb_download_data}/ukb1404_cal_chr1_v2_CURRENT.fixCol6.fam | sort > ${ukbb_sample_names_dir}/sort.all_samples.txt

echo -e "\t Making non-British European keep list\n"
comm -23 ${ukbb_sample_names_dir}/sort.all_samples.txt ${ukbb_sample_names_dir}/remove.non_british_european.txt > ${ukbb_sample_names_dir}/keep.non_british_european_122K.txt
