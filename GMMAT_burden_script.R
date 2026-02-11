args = commandArgs(trailingOnly=TRUE)

library(GMMAT)
library(SeqArray)

input_vcf <- args[1]
input_group_file <- args[2]
input_null_mod <- args[3]

output_gds <- sub("\\.vcf\\.gz$", ".gds", input_vcf)
output_file <- sub("\\.vcf\\.gz$", "_GMMAT-output.txt", input_vcf)

seqVCF2GDS(input_vcf, output_gds)

null_model <- readRDS(input_null_mod)

smmat_burden <- SMMAT(null_model, group.file = input_group_file, geno.file = output_gds, method = "davies", tests = "O")

output_table <- matrix(nrow = 1, ncol = 8)
output_table[1,1:4] <- args[1:4]
output_table[1,5:8] <- c(smmat_burden$n.variants, smmat_burden$B.pval, smmat_burden$S.pval, smmat_burden$O.pval)


write.table(output_table, file = output_file, quote=FALSE, sep="\t", row.names=FALSE, col.names = FALSE, na="")