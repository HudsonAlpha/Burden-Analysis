#!/bin/bash

#SBATCH -c 1
#SBATCH --mem=16G

if [ $# -ne 6 ] && [ $# -ne 7 ];
then
    echo "Usage: <filename> <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> [<arg7>]"
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
module load bcftools
module load cluster/java

# set variables for inputs
input_vcf_basename=$2
processing_dir_path=$3
filter_regions=$4
fam_input=$5
cov_input=$6

if [ $# -eq 7 ]; then
    kin_matrix=$7
else
    kin_matrix=""
fi

# path for skat R script
skat_script="/cluster/home/jtaylor/scripts/Burden_Analysis/skat_script.R"
skat_script_kin="/cluster/home/jtaylor/scripts/Burden_Analysis/skat_script_kinship.R"

scores="merged cadd fathmm genocanyon regulome eigen"

# make a directory for the gene
mkdir ${gene}
cd ${gene}

# set working dir
working_dir=$(pwd)

source /cluster/home/jtaylor/micromamba/etc/profile.d/micromamba.sh
micromamba activate /cluster/home/jtaylor/micromamba/envs/tools

slurm_mem_gb=$(echo "${SLURM_MEM_PER_NODE}/1024" | bc)

micromamba deactivate

# retrieve links and loops for the gene

awk -v gene="$gene" '($4 == gene || $6 == gene)' ${filter_regions} > ${gene}.bed

grep -w "${gene}" ${filter_regions} > ${gene}.bed

# find the chromosome for the gene
#chr=$(head -1 ${gene}.bed | awk '{ print $1 }')
chr=$(awk '{print $1}' ${gene}.bed | sort | uniq -c | sort -nr | head -1 | awk '{print $2}')

afs=(0.1 0.5 0.01 0.05 0.001)
cadd_scores=(10 20)

# 
for af in "${afs[@]}"; do
	for c_score in "${cadd_scores[@]}"; do 
	
		bcftools view -R ${gene}.bed ${processing_dir_path}/${chr}/${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}.vcf.gz -Ov \
			-o ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}.vcf --threads ${SLURM_JOB_CPUS_PER_NODE}
		
		bcftools view -R ${gene}.bed ${processing_dir_path}/${chr}/${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding.vcf.gz -Ov \
			-o ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}.vcf --threads ${SLURM_JOB_CPUS_PER_NODE}
		
		bcftools view -R ${gene}.bed ${processing_dir_path}/${chr}/${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF.vcf.gz -Ov \
			-o ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}.vcf --threads ${SLURM_JOB_CPUS_PER_NODE}
		
		
		
		bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}.vcf > \
			${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}.txt
		
		bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}.vcf > \
			${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}.txt
		
		bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}.vcf > \
			${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}.txt
		
		
		
		cat ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}.vcf | java -Xmx${slurm_mem_gb}G -jar /cluster/home/ncochran/bin/snpEff-4.3s/snpEff/SnpSift.jar extractFields - GEN[*].GT > \
			${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}_Gen.txt 
		sed -i 's/|/\//g; s/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s/1\/0/1/g; s/\.\/\./0/g; s/0\/\./0/g; s/\.\/0/0/g' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}_Gen.txt
		
		cat ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}.vcf | java -Xmx${slurm_mem_gb}G -jar /cluster/home/ncochran/bin/snpEff-4.3s/snpEff/SnpSift.jar extractFields - GEN[*].GT > \
			${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}_Gen.txt 
		sed -i 's/|/\//g; s/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s/1\/0/1/g; s/\.\/\./0/g; s/0\/\./0/g; s/\.\/0/0/g' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}_Gen.txt
		
		cat ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}.vcf | java -Xmx${slurm_mem_gb}G -jar /cluster/home/ncochran/bin/snpEff-4.3s/snpEff/SnpSift.jar extractFields - GEN[*].GT > \
			${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}_Gen.txt 
		sed -i 's/|/\//g; s/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s/1\/0/1/g; s/\.\/\./0/g; s/0\/\./0/g; s/\.\/0/0/g' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}_Gen.txt
		
		
		# add env for bc
		source /cluster/home/jtaylor/micromamba/etc/profile.d/micromamba.sh
		micromamba activate /cluster/home/jtaylor/micromamba/envs/tools
		
		
		awk '(NR>1)' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}_Gen.txt | /cluster/home/ncochran/bin/datamash-install/bin/datamash transpose | \
			sed 's/\t/+/g' | bc > ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}_wGen_xPose-Sums.txt
		
		awk '(NR>1)' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}_Gen.txt | /cluster/home/ncochran/bin/datamash-install/bin/datamash transpose | \
			sed 's/\t/+/g' | bc > ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}_wGen_xPose-Sums.txt
		
		awk '(NR>1)' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}_Gen.txt | /cluster/home/ncochran/bin/datamash-install/bin/datamash transpose | \
			sed 's/\t/+/g' | bc > ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}_wGen_xPose-Sums.txt
		
		
		micromamba deactivate
		
		
		# check file generated from last step. If 0, print 0; otherwise print 1, this ensures that individuals with more than one qualifying variant are only counted once
		awk '{if ($1=="0") print "0"; else print "1";}' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}_wGen_xPose-Sums.txt > \
			${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}_wGen_xPose-Sums-Collapse.txt
		
		awk '{if ($1=="0") print "0"; else print "1";}' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}_wGen_xPose-Sums.txt > \
			${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}_wGen_xPose-Sums-Collapse.txt
		
		awk '{if ($1=="0") print "0"; else print "1";}' ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}_wGen_xPose-Sums.txt > \
			${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}_wGen_xPose-Sums-Collapse.txt
	
		
		
		micromamba activate /cluster/home/jtaylor/micromamba/envs/skat
		
		if [ -n "${kin_matrix}" ]; then
			
			# run the skat R script
			Rscript ${skat_script_kin} ${fam_input} ${cov_input} ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}_wGen_xPose-Sums-Collapse.txt ${working_dir}/${gene}AF-${af}_CADD-${c_score}_all ${gene} ${kin_matrix}
		
			Rscript ${skat_script_kin} ${fam_input} ${cov_input} ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}_wGen_xPose-Sums-Collapse.txt ${working_dir}/${gene}_AF-${af}_CADD-${c_score}_coding ${gene} ${kin_matrix}
		
			Rscript ${skat_script_kin} ${fam_input} ${cov_input} ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}_wGen_xPose-Sums-Collapse.txt ${working_dir}/${gene}_AF-${af}_CADD-${c_score}_LOF ${gene} ${kin_matrix}
		
		else
			
			# run the skat R script
			Rscript ${skat_script} ${fam_input} ${cov_input} ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_${gene}_wGen_xPose-Sums-Collapse.txt ${working_dir}/${gene}AF-${af}_CADD-${c_score}_all ${gene}
		
			Rscript ${skat_script} ${fam_input} ${cov_input} ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_coding_${gene}_wGen_xPose-Sums-Collapse.txt ${working_dir}/${gene}_AF-${af}_CADD-${c_score}_coding ${gene}
		
			Rscript ${skat_script} ${fam_input} ${cov_input} ${input_vcf_basename}_${chr}_annotated_AF-${af}_CADD-${c_score}_LOF_${gene}_wGen_xPose-Sums-Collapse.txt ${working_dir}/${gene}_AF-${af}_CADD-${c_score}_LOF ${gene}
		
		fi
		
		micromamba deactivate
		
	done
	
done

# clean up temp files
rm *.vcf
rm *_Gen.txt
rm *_wGen_xPose-Sums.txt
