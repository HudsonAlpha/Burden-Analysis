#!/bin/bash

#SBATCH -c 1
#SBATCH --mem=8G

# check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <output dir> <processing dir"
    exit 1
fi

output_dir=$1
processing_dir=$2

:<<"COMMENT"
# create header for output file
echo -e "Fam file\tCovar set\tTest file\tGene\tSKAT p Val\tTotal vars\tCases with\tCases without\tControls With\tControls Without\tFET p\tOdds Ratio\tOR lower 95%\tOR upper 95%" > ${output_dir}/merged_burden_output.tsv

# combine all results to the output file
for file in ${processing_dir}/*/*_SKAT-output.txt; do
	cat ${file} >> ${output_dir}/merged_burden_output.tsv
done
COMMENT

# create header for output file
echo -e "test_vcf\tvariant_file\tnull_model\tgene\tnumber_variants\tburden_pval\tskat_pval\tskato_pval" > ${output_dir}/merged_burden_output.tsv

# combine all results to the output file
for file in ${processing_dir}/*/*_GMMAT-output.txt; do
	cat ${file} >> ${output_dir}/merged_burden_output.tsv
done