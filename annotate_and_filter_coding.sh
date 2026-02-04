#!/bin/bash

#SBATCH -p normal
#SBATCH -c 32
#SBATCH --mem=128G

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

module load bcftools
module load htslib
module load cluster/java
module load cluster/python/3.11.1
module load cluster/singularity

input_vcf_fn=$(basename ${input_vcf})
input_vcf_basename=${input_vcf_fn%.vcf.gz}

source /cluster/home/jtaylor/micromamba/etc/profile.d/micromamba.sh
micromamba activate /cluster/home/jtaylor/micromamba/envs/tools

slurm_mem_gb=$(echo "${SLURM_MEM_PER_NODE}/1024" | bc)

micromamba deactivate

static_config="/cluster/home/jtaylor/software/wgsa_v0.95/wgsa_static_config.txt"
wgsa_anno_head="/cluster/home/jtaylor/software/wgsa_v0.95/vcf_anno_header.txt"
wgsa_run_config="${input_vcf_basename}_wgsa_config.setting"
processing_script=" /cluster/home/jtaylor/scripts/Burden_Analysis/process_wgsa_output.py"
snp_sift="/cluster/home/ncochran/bin/snpEff_5.0/SnpSift.jar"

export CADD="/cluster/home/jtaylor/software/CADD-scripts"

working_dir=$(pwd)

# filter the VCF for the specified chromosome
bcftools view -R chr${chr}.bed ${input_vcf} -Oz -o chr${chr}/${input_vcf_basename}_filtered.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}

cd chr${chr}

tabix -p vcf ${input_vcf_basename}_filtered.vcf.gz

# decompose and normalize VCF with vt
#/cluster/software/vt-0.5772/bin/vt decompose ${input_vcf_basename}_filtered.vcf.gz -s -o + | \
#	/cluster/software/vt-0.5772/bin/vt normalize + -n -r /cluster/home/jtaylor/reference_files/hg38_asm5_alt/hg38.fa -o + | \
#	/cluster/software/vt-0.5772/bin/vt uniq + -o + | \
#	/cluster/software/vt-0.5772/bin/vt sort + | \
#	sed -e "s/chr//" > ${input_vcf_basename}_vt-noCHR.vcf

/cluster/software/vt-0.5772/bin/vt decompose ${input_vcf_basename}_filtered.vcf.gz -s -o temp1.vcf

/cluster/software/vt-0.5772/bin/vt normalize temp1.vcf -n -r /cluster/home/jtaylor/reference_files/hg38_asm5_alt/hg38.fa -o temp2.vcf

rm temp1.vcf

/cluster/software/vt-0.5772/bin/vt uniq temp2.vcf -o temp3.vcf

rm temp2.vcf

/cluster/software/vt-0.5772/bin/vt sort temp3.vcf -o temp4.vcf 

rm temp3.vcf

sed -e "s/chr//" temp4.vcf > ${input_vcf_basename}_vt-noCHR.vcf

rm temp4.vcf

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

: <<'END_COMMENT'
# move to CADD dir because CADD...
cd ${CADD}

# activate micromamba environment for CADD
micromamba activate /cluster/home/jtaylor/micromamba/envs/snakemake

# score VCF with CADD
#/cluster/home/ncochran/bin/CADD-scripts-master/CADD.sh -g GRCh38 -c ${SLURM_JOB_CPUS_PER_NODE} \
#	-o ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.tsv.gz ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.vcf.gz

# CADD 1.7 (no anno, CADD raw score and phred only)
#snakemake ${working_dir}/chr${chr}/${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.tsv.gz --sdm conda apptainer --apptainer-prefix /cluster/home/jtaylor/software/CADD-scripts/envs/apptainer \
#	--singularity-args "--bind /cluster:/cluster " --conda-prefix /cluster/home/jtaylor/software/CADD-scripts/envs/conda --cores 8 \
#	--configfile /cluster/home/jtaylor/software/CADD-scripts/config/config_GRCh38_v1.7_noanno.yml --snakefile /cluster/home/jtaylor/software/CADD-scripts/Snakefile -q


# CADD 1.7.2
./CADD.sh -g GRCh38 -c ${SLURM_JOB_CPUS_PER_NODE} -d -o ${working_dir}/chr${chr}/${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.tsv.gz ${working_dir}/chr${chr}/${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.tsv.gz


# deactivate environment
micromamba deactivate

# move back to chr directory
cd ${working_dir}/chr${chr}

tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.tsv.gz

# add 'chr' prefix and annotate VCF with CADD scores
bcftools annotate -a ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.tsv.gz -c CHROM,POS,REF,ALT,CADD_1.7_raw,CADD_1.7_phred \
	-h /cluster/home/jtaylor/reference_files/burden_analysis/cadd_head.txt -Ov ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.vcf.gz | sed '/^#/! s/^/chr/' \
	> ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7.vcf

END_COMMENT

# add 'chr' prefix and annotate VCF with CADD scores
bcftools annotate -a /cluster/home/jtaylor/software/CADD-scripts/data/prescored/GRCh38_v1.7/no_anno/whole_genome_SNVs.tsv.gz -c CHROM,POS,REF,ALT,CADD_1.7_raw,CADD_1.7_phred \
	-h /cluster/home/jtaylor/reference_files/burden_analysis/cadd_head.txt -Ov ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn.vcf.gz | sed '/^#/! s/^/chr/' \
	> ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7.vcf

bgzip -@${SLURM_JOB_CPUS_PER_NODE} ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7.vcf

tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7.vcf.gz

# convert missing genotypes to reference alleles
bcftools +missing2ref --threads ${SLURM_JOB_CPUS_PER_NODE} -Oz \
	-o ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7_m2r.vcf.gz ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7.vcf.gz

# index VCF after converting missing genotypes
tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7_m2r.vcf.gz

# fill VCF tags for allele counts, number, and frequency
bcftools +fill-tags ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7_m2r.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE} -Oz \
	-o ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7_m2r_FillTags.vcf.gz -- -t AN,AC,AF

tabix -p vcf ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7_m2r_FillTags.vcf.gz

# AFs and CADD scores to test
afs=(0.1 0.5 0.01 0.05 0.001)
cadd_scores=(10 20)

# filter vcf for AFs 
for af in "${afs[@]}"; do
	
	bcftools filter -i "INFO/Bravo_AF < ${af}" ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.7_m2r_FillTags.vcf.gz -Oz \
		-o ${input_vcf_basename}_chr${chr}_annotated_AF-${af}.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
	
	tabix -p vcf ${input_vcf_basename}_chr${chr}_annotated_AF-${af}.vcf.gz
	
	for c_score in "${cadd_scores[@]}"; do 
		
		bcftools filter -i "INFO/CADD_1.7_phred > ${c_score}" ${input_vcf_basename}_chr${chr}_annotated_AF-${af}.vcf.gz -Oz \
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
