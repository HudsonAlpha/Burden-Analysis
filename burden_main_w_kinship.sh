#!/bin/bash

#SBATCH -c 16
#SBATCH --mem=64G

# Function to display usage information
usage() {
    echo "Usage: $0 -v <input.vcf.gz> -f <input FAM> -c <input Cov> -b <coding or noncoding> [-k] [-h]"
    echo "Options:"
    echo "  -v  Input VCF file (gzipped)"
    echo "  -f  Input FAM file"
    echo "  -c  Input Covariates file"
    echo "  -b  Burden type ('coding' or 'noncoding')"
    echo "  -k  Use kinship matrix (flag, no input needed)"
		echo "  -n  Input Null Model rds"
    echo "  -h  Display this help message"
    exit 1
}

# Initialize kinship matrix flag as false
use_kin_matrix=false

# Parse command-line options using getopts
while getopts "v:f:c:b:k:h" opt; do
  case $opt in
    v) input_vcf="$OPTARG" ;;
    f) fam_input="$OPTARG" ;;
    c) cov_input="$OPTARG" ;;
    b) burden_type="$OPTARG" ;;
    k) use_kin_matrix=true
       kin_input="$OPTARG"
       ;;
		n) null_mod="$OPTARG" ;;
    h) usage ;;
    \?) usage ;;
  esac
done

# Check if required arguments are provided
if [ -z "$input_vcf" ] || [ -z "$fam_input" ] || [ -z "$cov_input" ] || [ -z "$burden_type" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Load modules
module load bcftools
module load htslib

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

: <<'END_COMMENT'
if [ "$use_kin_matrix" = true ]; then
	
	mkdir kinship
	kinship_matrix_path=$(realpath kinship)
	
	# make plink dataset for kinship matrix generation
	plink --vcf ${input_vcf} --const-fid --biallelic-only strict --make-bed --out ${input_vcf_basename} --memory ${SLURM_MEM_PER_NODE}
	
	awk '$6 != "NA"' ${fam_input} > temp.fam
	cat temp.fam | awk '{ print $1, $2 }' > temp_plink.filter
	
	plink --bfile ${input_vcf_basename} --make-bed --keep temp_plink.filter --out ${input_vcf_basename}_filtered --memory ${SLURM_MEM_PER_NODE}
	
	mv ${input_vcf_basename}_filtered.fam ${input_vcf_basename}_filtered_old.fam
	
	mv temp.fam ${input_vcf_basename}_filtered.fam
	
	/cluster/home/ncochran/bin/gemma-0.98.5-linux-static-AMD64 -bfile ${input_vcf_basename}_filtered -gk 2 -outdir "${kinship_matrix_path}/" -o ${input_vcf_basename}_filtered
	
	plink --bfile ${input_vcf_basename}_filtered --recode vcf-iid -out ${input_vcf_basename}_filtered --memory ${SLURM_MEM_PER_NODE}
	
	sed -i '/^#/! s/^/chr/' ${input_vcf_basename}_filtered.vcf
	
	bgzip -@ ${SLURM_JOB_CPUS_PER_NODE} ${input_vcf_basename}_filtered.vcf
	
	tabix -p vcf ${input_vcf_basename}_filtered.vcf.gz
	
	head -1 ${cov_input} > temp_cov_header.txt
	
	cat ${input_vcf_basename}_filtered.fam | awk '{ print $2 }' | while read line; do grep "${line}" ${cov_input}; done > temp_cov_body.txt
	
	cat temp_cov_header.txt temp_cov_body.txt > ${input_vcf_basename}_filtered.cov
	
	test_vcf=$(realpath ${input_vcf_basename}_filtered.vcf.gz)
	test_fam=$(realpath ${input_vcf_basename}_filtered.fam)
	test_cov=$(realpath ${input_vcf_basename}_filtered.cov)
	test_kin_matrix=$(realpath ${kinship_matrix_path}/${input_vcf_basename}_filtered.sXX.txt)
	
else
	
	test_vcf=${input_vcf}
	test_fam=${fam_input}
	test_cov=${cov_input}
	
fi

END_COMMENT


test_vcf=${input_vcf}
test_fam=${fam_input}
test_cov=${cov_input}
test_kin_matrix=${kin_input}

if [ $burden_type == "coding" ]; then
	
	filter_regions="/cluster/home/jtaylor/reference_files/burden_analysis/coding_regions_filt.bed"
	
	awk '{print $0 > $1".bed"}' ${filter_regions}
	
	# run all annotate and filter jobs for each chromosome
	af_job_id=$(sbatch --parsable --array [1-24]%6 -N 1 -n 1 --job-name "annotate_and_filter" \
		-o ${log_dir_path}/annotate_and_filter-%j_%A.out -e ${log_dir_path}/annotate_and_filter-%j_%A.err \
		/cluster/home/jtaylor/scripts/Burden_Analysis/annotate_and_filter_coding.sh /cluster/home/jtaylor/reference_files/burden_analysis/chr_file.txt ${test_vcf} ${log_dir_path})
	
	if [ "$use_kin_matrix" = true ]; then
		
		# run first batch of gene burden jobs
		gene_burden_1_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_1 \
			-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden_coding.sh \
			/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_1.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${test_fam} ${test_cov} ${null_mod} ${test_kin_matrix})

		# run second batch of gene burden jobs
		gene_burden_2_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_2 \
			-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden_coding.sh \
			/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_2.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${test_fam} ${test_cov} ${null_mod} ${test_kin_matrix})
		
	else	

		# run first batch of gene burden jobs
		gene_burden_1_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_1 \
			-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden_coding.sh \
			/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_1.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${test_fam} ${test_cov})

		# run second batch of gene burden jobs
		gene_burden_2_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_2 \
			-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden_coding.sh \
			/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_2.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${test_fam} ${test_cov})
				
	fi		
	
elif [ $burden_type == "noncoding" ]; then 
	
	filter_regions="/cluster/home/jtaylor/reference_files/burden_analysis/links_and_loops_sorted.bed"
	
	awk '{print $0 > $1".bed"}' ${filter_regions}
	
	# run all annotate and filter jobs for each chromosome
	af_job_id=$(sbatch --parsable --array [1-24]%6 -N 1 -n 1 --job-name "annotate_and_filter" \
		-o ${log_dir_path}/annotate_and_filter-%j_%A.out -e ${log_dir_path}/annotate_and_filter-%j_%A.err \
		/cluster/home/jtaylor/scripts/Burden_Analysis/annotate_and_filter.sh /cluster/home/jtaylor/reference_files/burden_analysis/chr_file.txt ${test_vcf} ${log_dir_path})

	# run first batch of gene burden jobs
	gene_burden_1_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_1 \
		-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden.sh \
		/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_1.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${test_fam} ${test_cov})

	# run second batch of gene burden jobs
	gene_burden_2_job=$(sbatch --parsable --array [1-10187]%20 -p normal -N 1 -n 1 --dependency=afterany:${af_job_id} --job-name gene_burden_2 \
		-o ${log_dir_path}/gene_burden-%j_%A.out -e ${log_dir_path}/gene_burden-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/gene_burden.sh \
		/cluster/home/jtaylor/reference_files/burden_analysis/gene_list_2.txt ${input_vcf_basename} ${processing_dir_path} ${filter_regions} ${test_fam} ${test_cov})

fi

# run job to combine gene burden output for all jobs		
sbatch -p normal --dependency=afterany:${gene_burden_1_job},${gene_burden_2_job} --job-name combine_burden_results -o ${log_dir_path}/combine_burden.out \
	-e ${log_dir_path}/combine_burden.err /cluster/home/jtaylor/scripts/Burden_Analysis/combine_burden_results.sh ${starting_dir_path} ${processing_dir_path}