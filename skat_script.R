args = commandArgs(trailingOnly=TRUE)

library(SKAT)

# expects argument 1 to be a fam file, argument 2 to be a covariate file in plink format, and argument 3 to be a burden set
# the burden set as well as the phenotype designation should be binary and 0/1 coded
FAM_Cov<-Read_Plink_FAM_Cov(args[1],args[2],Is.binary=TRUE,flag1=1)
BurdenSet <- as.matrix(read.table(args[3], header=FALSE, sep ="\t"))

Pheno<-FAM_Cov$Phenotype

# determine the covariates' column names, excluding the phenotype column if it's included
covariate_columns <- setdiff(names(FAM_Cov), c("FID", "IID", "PID", "MID", "Phenotype"))

# Remove covariate columns with all values being NA
#covariate_columns <- covariate_columns[sapply(FAM_Cov[covariate_columns], function(x) !all(is.na(x)))]

# dynamically create variables for each covariate
for (covariate in covariate_columns) {
  assign(covariate, FAM_Cov[[covariate]])
}

# create the model formula string dynamically
covariates_formula_part <- paste(covariate_columns, collapse=" + ")
formula_str <- paste("Pheno ~", covariates_formula_part)

cat("Null model formula:", formula_str, "\n")

# convert the formula string into a formula object
formula <- as.formula(formula_str)

# fit the SKAT null model using the dynamically created formula
obj <- SKAT_Null_Model(formula, data = FAM_Cov, out_type = "D", Adjustment = FALSE)

# perform SKAT with robust method
SKATout<-SKATBinary_Robust(BurdenSet,obj,method="Burden")
SKATout$p.value
SKATout$mac
SKATout$param$n.marker


# create a matrix to store results
BigTable <- matrix(nrow = 1, ncol = 7)
BigTable[1,1:3] <- args[1:3]
BigTable[1,4] <- args[5]
BigTable[1,5:7] <- c(SKATout$p.value, SKATout$mac, SKATout$param$n.marker)

# SKAT part is done, now get effect size and run Fisher's Exact Test
TotalControls <- as.data.frame(table(Pheno))[1,2]
TotalCases <- as.data.frame(table(Pheno))[2,2]
PhenoControls <- (Pheno-1)*-1
CasesWithVar <- as.data.frame(table(Pheno*BurdenSet))[2,2]
CasesWithVar[is.na(CasesWithVar)] <- 0
CasesWithoutVar <- TotalCases-CasesWithVar
ControlsWithVar <- as.data.frame(table(PhenoControls*BurdenSet))[2,2]
ControlsWithVar[is.na(ControlsWithVar)] <- 0
ControlsWithoutVar <- TotalControls-ControlsWithVar

# function to run Fisher's Exact Test on a given row of data
row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           p_val = f$p.value,
           or = f$estimate[[1]],
           or_ll = f$conf.int[1],
           or_ul = f$conf.int[2]))
}

# prepare data for Fisher's Exact Test
TestData <- matrix(nrow = 1, ncol = 4)
TestData[1,1:4] = c(CasesWithVar, CasesWithoutVar, ControlsWithVar, ControlsWithoutVar)
colnames(TestData) <- c('CaW', 'CaWo', 'CoW', 'CoWo')

# apply Fisher's test to data
FET <- t(apply(TestData, 1, row_fisher))

# combine all results
OutputTable <- cbind(BigTable,FET)
OutputName <- paste(args[4],"_SKAT-output.txt",sep = "")

write.table(OutputTable,file=OutputName,quote=FALSE,sep="\t",row.names=FALSE,col.names = FALSE,na="")