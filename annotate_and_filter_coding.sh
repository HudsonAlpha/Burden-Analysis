#!/bin/bash

#SBATCH -c 8
#SBATCH --mem=64G

if [ $# -ne 3 ];
then
    echo "Usage: <chr file> <input vcf> <log dir>"
    exit 1
fi

if [ ! -f $1 ];
then
    echo "$1 file does not exist"
    exit 1
fi

chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p $1)
echo $chr

input_vcf=$2
log_dir=$3

module load cluster/bcftools
module load cluster/htslib
module load cluster/java
module load cluster/python/3.7.6

input_vcf_fn=$(basename ${input_vcf})
input_vcf_basename=${input_vcf_fn%.vcf.gz}

slurm_mem_gb=$(echo "${SLURM_MEM_PER_NODE}/1024" | bc)

static_config="/cluster/home/jtaylor/software/wgsa_v0.95/wgsa_static_config.txt"
wgsa_anno_head="/cluster/home/jtaylor/software/wgsa_v0.95/vcf_anno_header.txt"
wgsa_run_config="${input_vcf_basename}_wgsa_config.setting"
processing_script=" /cluster/home/jtaylor/scripts/Burden_Analysis/process_wgsa_output.py"
snp_sift="/cluster/home/ncochran/bin/snpEff_5.0/SnpSift.jar"

working_dir=$(pwd)

# filter the VCF for the specified chromosome
bcftools view -R chr${chr}.bed ${input_vcf} -Oz -o chr${chr}/${input_vcf_basename}_filtered.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}

cd chr${chr}

tabix -p vcf ${input_vcf_basename}_filtered.vcf.gz

# decompose and normalize VCF with vt
/cluster/software/vt-0.5772/bin/vt decompose ${input_vcf_basename}_filtered.vcf.gz -s -o + | \
	/cluster/software/vt-0.5772/bin/vt normalize + -r /cluster/home/ncochran/Scripts/hg38.fa -o + | \
	/cluster/software/vt-0.5772/bin/vt uniq + -o + | \
	/cluster/software/vt-0.5772/bin/vt sort + | \
	sed -e "s/chr//" > ${input_vcf_basename}_vt-noCHR.vcf

# index decomposed and normalized VCF
bgzip -@${SLURM_JOB_CPUS_PER_NODE} ${input_vcf_basename}_vt-noCHR.vcf
tabix -p vcf ${input_vcf_basename}_vt-noCHR.vcf.gz

# annotate VCF with dbSNP IDs
bcftools annotate -a /cluster/home/jtaylor/reference_files/dbSNP/dbSNP_156.vcf.gz -c ID --threads ${SLURM_JOB_CPUS_PER_NODE} \
	-Oz -o ${input_vcf_basename}_dbSNP-156.vcf.gz ${input_vcf_basename}_vt-noCHR.vcf.gz

tabix -p vcf ${input_vcf_basename}_dbSNP-156.vcf.gz

# annotate VCF with TOPMed info
bcftools annotate -a /cluster/home/jtaylor/reference_files/burden_analysis/chrALL.BRAVO_TOPMed_Freeze_8_NOchr_Renamed.vcf.gz \
	-c INFO --threads ${SLURM_JOB_CPUS_PER_NODE} -Oz -o ${input_vcf_basename}_dbSNP-156_wTOPMed.vcf.gz ${input_vcf_basename}_dbSNP-156.vcf.gz

# index VCF annotated with TOPMed info
tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed.vcf.gz

# annotate VCF with snpEff predictions
java -Xmx64G -jar /cluster/home/ncochran/bin/snpEff_5.0/snpEff.jar -s ${input_vcf_basename}_snpeff_stats.html \
	GRCh38.99 ${input_vcf_basename}_dbSNP-156_wTOPMed.vcf.gz > ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.vcf
	
bgzip -@${SLURM_JOB_CPUS_PER_NODE} ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.vcf

# index VCF annotated with snpEff predictions
tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.vcf.gz

# activate micromamba environment for CADD
source /cluster/home/jtaylor/micromamba/etc/profile.d/micromamba.sh
micromamba activate /cluster/home/jtaylor/micromamba/envs/snakemake

# score VCF with CADD
/cluster/home/ncochran/bin/CADD-scripts-master/CADD.sh -g GRCh38 -c ${SLURM_JOB_CPUS_PER_NODE} \
	-o ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.tsv.gz ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.vcf.gz

# deactivate environment (check if micromamba or conda)
micromamba deactivate

tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.tsv.gz

# add 'chr' prefix and annotate VCF with CADD scores
bcftools annotate -a ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.tsv.gz -c CHROM,POS,REF,ALT,CADD_1.6_raw,CADD_1.6_phred \
	-h /cluster/home/ncochran/Scripts/cadd_head.txt -Ov ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.vcf.gz | sed '/^#/! s/^/chr/' \
	> ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6.vcf

bgzip -@${SLURM_JOB_CPUS_PER_NODE} ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6.vcf

tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6.vcf.gz

# convert missing genotypes to reference alleles
bcftools +missing2ref --threads ${SLURM_JOB_CPUS_PER_NODE} -Oz \
	-o ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6_m2r.vcf.gz ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6.vcf.gz

# index VCF after converting missing genotypes
tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6_m2r.vcf.gz

# fill VCF tags for allele counts, number, and frequency
bcftools +fill-tags ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6_m2r.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE} -Oz \
	-o ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6_m2r_FillTags.vcf.gz -- -t AN,AC,AF

tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6_m2r_FillTags.vcf.gz

# AFs and CADD scores to test
afs=(0.01 0.001 0.0001)
cadd_scores=(10 20)

# filter vcf for AFs 
for af in "${afs[@]}"; do
	
	bcftools filter -i "INFO/Bravo_AF < ${af}" ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6_m2r_FillTags.vcf.gz -Oz \
		-o ${input_vcf_basename}_chr${chr}_annotated_AF-${af}.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
	
	tabix -p vcf ${input_vcf_basename}_chr${chr}_annotated_AF-${af}.vcf.gz
	
	for c_score in "${cadd_scores[@]}"; do 
		
		bcftools filter -i "INFO/CADD_1.6_phred > ${c_score}" ${input_vcf_basename}_chr${chr}_annotated_AF-${af}.vcf.gz -Oz \
			-o ${input_vcf_basename}_chr${chr}_annotated_AF-${af}_CADD-${c_score}.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
		
		tabix -p vcf ${input_vcf_basename}_chr${chr}_annotated_AF-${af}_CADD-${c_score}.vcf.gz
		
		# filter vcf for protein coding variants only
		java -Xmx${slurm_mem_gb}G -jar ${snp_sift} filter "((ANN[0].IMPACT = 'HIGH') | (ANN[0].IMPACT = 'MODERATE'))" ${input_vcf_basename}_chr${chr}_annotated_AF-${af}_CADD-${c_score}.vcf.gz > \
			${input_vcf_basename}_chr${chr}_annotated_AF-${af}_CADD-${c_score}_coding.vcf
		
		bgzip -@${SLURM_JOB_CPUS_PER_NODE} ${input_vcf_basename}_chr${chr}_annotated_AF-${af}_CADD-${c_score}_coding.vcf
		
		tabix -p vcf ${input_vcf_basename}_chr${chr}_annotated_AF-${af}_CADD-${c_score}_coding.vcf.gz
		
		bcftools filter -i 'INFO/LOF ~ ".*"' ${input_vcf_basename}_chr${chr}_annotated_AF-${af}_CADD-${c_score}.vcf.gz -Oz \
			-o ${input_vcf_basename}_chr${chr}_annotated_AF-${af}_CADD-${c_score}_LOF.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
		
		tabix -p vcf ${input_vcf_basename}_chr${chr}_annotated_AF-${af}_CADD-${c_score}_LOF.vcf.gz
		
	done
	
done
