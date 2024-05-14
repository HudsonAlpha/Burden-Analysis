#!/bin/bash


***********************************
******* Still needs testing *******
***********************************


#SBATCH -c 1
#SBATCH --mem=16G

if [ $# -ne 6 ];
then
    echo "Usage: <filename>"
    exit 1
fi

if [ ! -f $1 ];
then
    echo "$1 file does not exist"
    exit 1
fi

# set gene for burden analysis
gene=$(sed -n ${SLURM_ARRAY_TASK_ID}p $1)
echo $gene

# load modules
module load cluster/bcftools
module load cluster/java

# set variables for inputs
input_vcf_basename=$2
processing_dir_path=$3
links_loops=$4
fam_input=$5
cov_input=$6

# path for skat R script
skat_script="/cluster/home/jtaylor/scripts/Burden_Analysis/skat_script.R"

scores="merged cadd fathmm genocanyon regulome eigen"

# make a directory for the gene
mkdir ${gene}
cd ${gene}

# set working dir
working_dir=$(pwd)

# retrieve links and loops for the gene
grep "${gene}" ${links_loops} > ${gene}.bed

# find the chromosome for the gene
chr=$(head -1 ${gene}.bed | awk '{ print $1 }')

# set path for vcf
merged_vcf="${processing_dir_path}/${chr}/${input_vcf_basename}_all_annotated_filtered_${chr}.vcf.gz"
cadd_vcf="${processing_dir_path}/${chr}/${input_vcf_basename}_noncoding_${chr}_annotated_filtered_AF-0.0001_impact_CADD-15.vcf.gz"
fathmm_vcf="${processing_dir_path}/${chr}/${input_vcf_basename}_noncoding_${chr}_annotated_filtered_AF-0.0001_impact_FXF-0.7.vcf.gz"
genocanyon_vcf="${processing_dir_path}/${chr}/${input_vcf_basename}_noncoding_${chr}_annotated_filtered_AF-0.0001_impact_GC-0.7.vcf.gz"
regulome_vcf="${processing_dir_path}/${chr}/${input_vcf_basename}_noncoding_${chr}_annotated_filtered_AF-0.0001_impact_RD.vcf.gz"
eigen_vcf="${processing_dir_path}/${chr}/${input_vcf_basename}_noncoding_${chr}_annotated_filtered_AF-0.0001_impact_Eigen-15.vcf.gz"

# subset vcf to snps that fall in associated regions for the gene
bcftools view -R ${gene}.bed ${merged_vcf} -Ov -o ${input_vcf_basename}_merged_${gene}.vcf --threads ${SLURM_JOB_CPUS_PER_NODE}
bcftools view -R ${gene}.bed ${cadd_vcf} -Ov -o ${input_vcf_basename}_cadd_${gene}.vcf --threads ${SLURM_JOB_CPUS_PER_NODE}
bcftools view -R ${gene}.bed ${fathmm_vcf} -Ov -o ${input_vcf_basename}_fathmm_${gene}.vcf --threads ${SLURM_JOB_CPUS_PER_NODE}
bcftools view -R ${gene}.bed ${genocanyon_vcf} -Ov -o ${input_vcf_basename}_genocanyon_${gene}.vcf --threads ${SLURM_JOB_CPUS_PER_NODE}
bcftools view -R ${gene}.bed ${regulome_vcf} -Ov -o ${input_vcf_basename}_regulome_${gene}.vcf --threads ${SLURM_JOB_CPUS_PER_NODE}
bcftools view -R ${gene}.bed ${eigen_vcf} -Ov -o ${input_vcf_basename}_eigen_${gene}.vcf --threads ${SLURM_JOB_CPUS_PER_NODE}

for score in ${scores}; do
	
	# extract genotype information from the VCF using SnpSift
	cat ${input_vcf_basename}_${score}_${gene}.vcf | java -Xmx16G -jar /cluster/home/ncochran/bin/snpEff-4.3s/snpEff/SnpSift.jar extractFields - GEN[*].GT > ${input_vcf_basename}_${score}_${gene}_Gen.txt 
	sed -i 's/|/\//g;s/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s/1\/0/1/g; s/\.\/\./0/g; s/0\/\./0/g; s/\.\/0/0/g' ${input_vcf_basename}_${score}_${gene}_Gen.txt
	
done

for score in ${scores}; do
	
	# transpose, replace tabs with +, and find the sum 
	awk '(NR>1)' ${input_vcf_basename}_${score}_${gene}_Gen.txt | /cluster/home/ncochran/bin/datamash-install/bin/datamash transpose | \
		sed 's/\t/+/g' | bc > ${input_vcf_basename}_${score}_${gene}_wGen_xPose-Sums.txt
	
done

for score in ${scores}; do
	
	# check file generated from last step. If 0, print 0; otherwise print 1, this ensures that individuals with more than one qualifying variant are only counted once
	awk '{if ($1=="0") print "0"; else print "1";}' ${input_vcf_basename}_${score}_${gene}_wGen_xPose-Sums.txt > ${input_vcf_basename}_${score}_${gene}_wGen_xPose-Sums-Collapse.txt
	
done

# clean up temp files
rm *.vcf
rm *_Gen.txt
rm *_wGen_xPose-Sums.txt

# activate env for skat
source /cluster/home/jtaylor/micromamba/etc/profile.d/micromamba.sh
micromamba activate /cluster/home/jtaylor/micromamba/envs/skat

for score in ${scores}; do
	
	# run the skat R script
	Rscript ${skat_script} ${fam_input} ${cov_input} ${input_vcf_basename}_${score}_${gene}_wGen_xPose-Sums-Collapse.txt ${working_dir}/${gene}_${score} ${gene}
	
done

micromamba deactivate