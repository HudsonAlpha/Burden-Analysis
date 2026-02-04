#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=64G

if [ $# -ne 2 ];
then
    echo "Usage: <filename>"
    exit 1
fi

if [ ! -f $1 ];
then
    echo "$1 file does not exist"
    exit 1
fi

filt_type=$(sed -n ${SLURM_ARRAY_TASK_ID}p $1)
echo $filt_type

module load cluster/bcftools
module load cluster/htslib

input_vcf=$2
input_vcf_basename=${input_vcf%.vcf.gz}

# filter input vcf for only variants with CADD phred > 15
if [[ $filt_type == "cadd" ]]
then
	
	bcftools filter -i 'INFO/CADD_1.6_phred > 15' ${input_vcf} -Oz -o ${input_vcf_basename}_CADD-15.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
	tabix -p vcf ${input_vcf_basename}_CADD-15.vcf.gz

# filter input vcf for only variants with genocanyon score > 0.5	
elif [[ $filt_type == "genocanyon" ]]
then
	
	bcftools filter -i 'INFO/GenoCanyon_score > 0.7' ${input_vcf} -Oz -o ${input_vcf_basename}_GC-0.7.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
	tabix -p vcf ${input_vcf_basename}_GC-0.7.vcf.gz

# filter input vcf for only varinats with a qualifying regulomedb score	
elif [[ $filt_type == "regulomedb" ]]
then
	
	filt1="INFO/RegulomeDB_score = \"1a\" || INFO/RegulomeDB_score = \"1b\" || INFO/RegulomeDB_score = \"1c\" || INFO/RegulomeDB_score = \"1d\""
	filt2="INFO/RegulomeDB_score = \"1e\" || INFO/RegulomeDB_score = \"1f\" || INFO/RegulomeDB_score = \"2a\" || INFO/RegulomeDB_score = \"2b\""
	filt3="INFO/RegulomeDB_score = \"2c\" || INFO/RegulomeDB_score = \"3a\" || INFO/RegulomeDB_score = \"3b\""

	bcftools filter -i "${filt1} || ${filt2} || ${filt3}" ${input_vcf} -Oz -o ${input_vcf_basename}_RD.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
	tabix -p vcf ${input_vcf_basename}_RD.vcf.gz

# filter input vcf for only variants with fathmmxf score > 0.5	
elif [[ $filt_type == "fathmmxf" ]]
then
	
	bcftools filter -i 'INFO/fathmmXF_score > 0.7' ${input_vcf} -Oz -o ${input_vcf_basename}_FXF-0.7.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
	tabix -p vcf ${input_vcf_basename}_FXF-0.7.vcf.gz

# filter input vcf for only variants with eigen phred > 15	
elif [[ $filt_type == "eigen" ]]
then
	
	bcftools filter -i 'INFO/Eigen_phred > 15' ${input_vcf} -Oz -o ${input_vcf_basename}_Eigen-15.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
	tabix -p vcf ${input_vcf_basename}_Eigen-15.vcf.gz
	
else
	
	echo "Not a valid filter type."
	
fi

