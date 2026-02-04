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
module load cluster/python/3.11.1

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
	/cluster/software/vt-0.5772/bin/vt normalize + -n -r /cluster/home/ncochran/Scripts/hg38.fa -o + | \
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

### wgsa

# make tmp and work dir for wgsa
mkdir tmp
mkdir work

# write config file for WGSA
echo "input file name: ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6_m2r_FillTags.vcf.gz" > temp_config.txt
echo "output file name: ${input_vcf_basename}_annotated" >> temp_config.txt
echo "resources dir: /cluster/home/jtaylor/software/wgsa_v0.95/resources/" >> temp_config.txt
echo "annovar dir: /cluster/home/jtaylor/software/wgsa/annovar20220328/annovar/" >> temp_config.txt
echo "snpeff dir: /cluster/home/jtaylor/software/wgsa/snpeff/snpEff/" >> temp_config.txt
echo "vep dir: /cluster/home/jtaylor/software/wgsa/vep/ensembl-vep-release-106/" >> temp_config.txt
echo ".vep dir: /cluster/home/jtaylor/software/wgsa/.vep/" >> temp_config.txt
echo "tmp dir: ${working_dir}/chr${chr}/tmp/" >> temp_config.txt
echo "work dir: ${working_dir}/chr${chr}/work/" >> temp_config.txt
cat temp_config.txt ${static_config} > ${wgsa_run_config}

# run initial step for WGSA
echo -e "Understand" | java -cp /cluster/home/jtaylor/software/wgsa_v0.95/ WGSA095 ${wgsa_run_config} \
	-m ${slurm_mem_gb} -t ${SLURM_JOB_CPUS_PER_NODE} -v hg38 -i vcf

# make WGSA config script executable
chmod +x ${input_vcf_basename}_wgsa_config.setting.sh

# run WGSA
bash ${input_vcf_basename}_wgsa_config.setting.sh

# process snps and indels from WGSA to create annotation file
python3 ${processing_script} ${input_vcf_basename}_annotated.snp.gz ${input_vcf_basename}_annotated.indel.gz ${input_vcf_basename}_annotated.wgsa.tsv

# zip the new annotation file
bgzip -@${SLURM_JOB_CPUS_PER_NODE} ${input_vcf_basename}_annotated.wgsa.tsv

# index annotation file
tabix -p vcf ${input_vcf_basename}_annotated.wgsa.tsv.gz

# annotate vcf with WGSA scores
bcftools annotate -a ${input_vcf_basename}_annotated.wgsa.tsv.gz \
	-c CHROM,POS,REF,ALT,GenoCanyon_score,RegulomeDB_score,fathmmXF_score,Eigen_phred \
	-h ${wgsa_anno_head} --threads ${SLURM_JOB_CPUS_PER_NODE} -Oz \
	-o ${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered.vcf.gz ${input_vcf_basename}_dbSNP-156_wTOPMed_wAnn_wCADD-1.6_m2r_FillTags.vcf.gz

# index vcf
tabix -p vcf ${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered.vcf.gz

# filter vcf for BRAVO af of 1:10,000
bcftools filter -i 'INFO/Bravo_AF < 0.0001' ${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered.vcf.gz -Oz \
	-o ${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}

#index vcf
tabix -p vcf ${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001.vcf.gz

# filter vcf based on impact to provide non-coding variant only
java -Xmx${slurm_mem_gb}G -jar ${snp_sift} filter "!((ANN[0].IMPACT = 'HIGH') | (ANN[0].IMPACT = 'MODERATE'))" ${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001.vcf.gz > \
	${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001_impact.vcf

# gzip new vcf
bgzip -@${SLURM_JOB_CPUS_PER_NODE} ${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001_impact.vcf

# index new vcf
tabix -p vcf ${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001_impact.vcf.gz

# submit job array of individual score filters
filt_job_id=$(sbatch --parsable --array [1-5]%3 -p normal -N 1 -n 1 --job-name "filter_annotations" -o ${log_dir}/filter_anno-%j_%A.out \
	-e ${log_dir}/filter_anno-%j_%A.err /cluster/home/jtaylor/scripts/Burden_Analysis/filter_annotations.sh \
	/cluster/home/jtaylor/reference_files/burden_analysis/filt_file.txt ${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001_impact.vcf.gz)

# submit job to wait for all filter jobs
wait_job_id=$(sbatch --parsable --job-name waiting_chr${chr} --dependency=afterany:${filt_job_id} --wrap="echo 'All jobs completed'")

# wait for filter jobs to finish
echo "Waiting for all jobs to complete..."
while squeue -j $wait_job_id | grep -q $wait_job_id; do
	sleep 100
done

echo "All SLURM jobs have completed."

# create a merge list of vcfs
echo "${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001_impact_CADD-15.vcf.gz" > merge.list
echo "${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001_impact_GC-0.7.vcf.gz" >> merge.list
echo "${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001_impact_RD.vcf.gz" >> merge.list
echo "${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001_impact_FXF-0.7.vcf.gz" >> merge.list
echo "${input_vcf_basename}_noncoding_chr${chr}_annotated_filtered_AF-0.0001_impact_Eigen-15.vcf.gz" >> merge.list

# concatenate all filter vcfs
bcftools concat -f merge.list -a -D -Oz -o ${input_vcf_basename}_all_annotated_filtered_chr${chr}.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}

# index vcf
tabix -p vcf ${input_vcf_basename}_all_annotated_filtered_chr${chr}.vcf.gz