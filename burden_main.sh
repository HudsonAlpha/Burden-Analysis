#!/bin/bash

#SBATCH -c 1
#SBATCH --mem=16G

# check if the correct number of arguments is provided
if [ "$#" -ne 4 ] && [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input.vcf.gz> <input FAM> <input Cov> <coding or noncoding> [<kinship matrix>]"
    exit 1
fi

# load modules
module load cluster/bcftools
module load cluster/htslib

# inputs
input_vcf=$1
fam_input=$2
cov_input=$3
burden_type=$4
kin_matrix=$5

# vcf basename
input_vcf_fn=$(basename ${input_vcf})
input_vcf_basename=${input_vcf_fn%.vcf.gz}

# chr variable
chromosomes=({1..22} X Y)

# create a directory for log files
log_dir="burden_logs_$(date +%Y%m%d)"
mkdir -p "$log_dir"
log_dir_path=$(realpath "$log_dir")

# set starting directory
starting_dir_path=$(pwd)

# create processing directory
mkdir -p processing
processing_dir_path=$(realpath processing)
cd ${processing_dir_path}

for chr in "${chromosomes[@]}"; do
	mkdir chr${chr}
done

# ensure the vcf is indexed
if [ ! -f "${input_vcf}.tbi" ]; then
	echo "Index file for '${input_vcf_fn}' not found. Indexing now..."

	tabix -p vcf ${input_vcf}

	if [ $? -eq 0 ]; then
		echo "VCF indexed."
	else
		echo "Indexing failed."
		exit 1
	fi
else
	echo "Index file for ${input_vcf} exists."
fi


if [ $burden_type == "coding" ]; then
	
	filter_regions="/cluster/home/jtaylor/reference_files/burden_analysis/coding_regions_filt.bed"
	
	awk '{print $0 > $1".bed"}' ${filter_regions}
	
	# run all annotate and filter jobs for each chromosome
	af_job_id=$(sbatch --parsable --array [1-24]%6 -p normal -N 1 -n 1 --job-name "annotate_and_filter" \
		-o ${log_dir_path}/annotate_and_filter-%j_%A.out -e ${log_dir_path}/annotate_and_filter-%j_%A.err \
		/cluster/home/jtaylor/scripts/Burden_Analysis/annotate_and_filter_coding.sh /cluster/home/jtaylor/reference_files/burden_analysis/chr_file.txt ${input_vcf} ${log_dir_path})
	
	if [ -n "${kin_matrix}" ]; then
		
		# run first batch of gene burden jobs
		gene_burden_1_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_1 \
			-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden_coding.sh \
			/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_1.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${fam_input} ${cov_input} ${kin_matrix})

		# run second batch of gene burden jobs
		gene_burden_2_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_2 \
			-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden_coding.sh \
			/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_2.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${fam_input} ${cov_input} ${kin_matrix})
		
	else	

		# run first batch of gene burden jobs
		gene_burden_1_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_1 \
			-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden_coding.sh \
			/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_1.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${fam_input} ${cov_input})

		# run second batch of gene burden jobs
		gene_burden_2_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_2 \
			-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden_coding.sh \
			/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_2.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${fam_input} ${cov_input})
				
	fi		
	
elif [ $burden_type == "noncoding" ]; then 
	
	filter_regions="/cluster/home/jtaylor/reference_files/burden_analysis/links_and_loops_sorted.bed"
	
	awk '{print $0 > $1".bed"}' ${filter_regions}
	
	# run all annotate and filter jobs for each chromosome
	af_job_id=$(sbatch --parsable --array [1-24]%6 -p normal -N 1 -n 1 --job-name "annotate_and_filter" \
		-o ${log_dir_path}/annotate_and_filter-%j_%A.out -e ${log_dir_path}/annotate_and_filter-%j_%A.err \
		/cluster/home/jtaylor/scripts/Burden_Analysis/annotate_and_filter.sh /cluster/home/jtaylor/reference_files/burden_analysis/chr_file.txt ${input_vcf} ${log_dir_path})

	# run first batch of gene burden jobs
	gene_burden_1_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_1 \
		-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden.sh \
		/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_1.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${fam_input} ${cov_input})

	# run second batch of gene burden jobs
	gene_burden_2_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_2 \
		-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden.sh \
		/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_2.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${fam_input} ${cov_input})

fi

# run job to combine gene burden output for all jobs		
sbatch -p normal --dependency=afterany:${gene_burden_1_job},${gene_burden_2_job} --job-name combine_burden_results -o ${log_dir_path}/combine_burden.out \
	-e ${log_dir_path}/combine_burden.err /cluster/home/jtaylor/scripts/Burden_Analysis/combine_burden_results.sh ${starting_dir_path} ${processing_dir_path}