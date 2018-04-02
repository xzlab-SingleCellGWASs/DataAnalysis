###This program aims to processing dataset, including phenotye and QC. 

library(data.table)
library(dplyr)
library(pipeR)

#####Input and output files 
###Input
##sqc data: ukb_sqc_v2.txt
##sqc head data: ukb_sqc_v2_head.txt (in the software)
##fam data: ukb30186_cal_chrXX_v2_s488363.fam (Any fam file is OK)
##phenotype: ukb10683.csv (RData)
##phenotype for analysis
###Output
##phenotype output: Selected samples of analyzing pheno
##sample output: Sample information(eid, famid1, famid2)  
##sqc output:

args <- commandArgs(TRUE)
sqcFile <- args[1]
hdrFile <- args[2]
famFile <- args[3]
phenoFile <- args[4]
nlyFile <- args[5]
phenoOut < args[6]
sampleOut <- args[7]
sqcOut < args[8]


#########STEP1: combine the sqc and fam
##sqc
##fam
##header of sqc
##Numbers of row of the sqc and fam should be the same.
print("STEP1: combine the sqc and fam...")
#File input
sqc <- fread(sqcFile, stringsAsFactors = F)
sqc <- data.frame(sqc)
header <- read.table(hdrFile, stringsAsFactors = F)
header <- data.frame(header)
fam <- fread(famFile, stringsAsFactors = F)
fam <- data.frame(fam)

#
if (dim(fam)[1] != dim(sqc)[1]){
 
  stop("ERROR: The number of rows is not the same!")
}else{
  
  sqc <- cbind.data.frame(fam[1], sqc)
  names(sqc) <- c("eid", header[, 1])
  cat("QC sample:", dim(sqc)[1], "\n")
}

#########STEP2: select the samples

##in.Phasing.Input.chr1_22==1
##in.white.British.ancestry.subset==1
##used.in.pca.calculation==1
##excess.relatives==0
##putative.sex.chromosome.aneuploidy==0
##eid > 0
print("STEP2: select the samples...")
##Process the QC results
#Genetyping
sqc.i <- sqc[which(sqc$in.Phasing.Input.chr1_22==1), ]
eid.i <- sqc.i[, 1]
cat("Genotyping success:", dim(sqc.i)[1], "\n")

#Others
cnd.i <- sqc.i$in.white.British.ancestry.subset == 1
cat("White British ancestry subset:", sum(cnd.i), "\n")
cnd.ii <- sqc.i$excess.relatives == 0 
cat("Excess relatives:", sum(!cnd.ii), "\n")
cnd.iii <- sqc.i$putative.sex.chromosome.aneuploidy == 0 
cat("Sex chromosome aneuploidy:", sum(!cnd.iii), "\n")
cnd.ix <- sqc.i$used.in.pca.calculation == 1
cat("Used in PCA calculation:", sum(cnd.ix), "\n")
cnd.x <- sqc.i$eid > 0
cat("Redacted:", sum(!cnd.x), "\n")

#Index of selected samples
cnd <- cnd.i & cnd.ii & cnd.iii & cnd.ix & cnd.x
cat("Samples Remaining:", sum(cnd), "\n")
write.csv(as.numeric(cnd), cndOut, row.names = F)


#########STEP3: phenotype data
print("STEP3: phenotype data...")
load(phenoFile)
label <- read.table(nlyFile, header = F)

#Select and sort the phenotype data
pheno <- ukb[which(ukb$eid %in% eid.i), ]
idx <- match(eid.i, pheno[, 1])
pheno <- pheno[idx, ]

if (all(label[, 1] %in% colnames(pheno)) == F){
  
  stop("ERROR: The labels of phenotype are wrong!")
}
pheno <- pheno[, c(1, which(colnames(pheno) %in% label[, 1]))]
pheno[which(is.na(pheno[, 1])), 1] <- c(-1:-14)

#
if (nrow(label) == 1){
  pheno[, 2] <- qqnorm(pheno[, 2], plot.it = F)$x
}else{
  pheno[, c(2: ncol(pheno))] <- apply(pheno[, c(2: ncol(pheno))], 2, function(a) 
                                      qqnorm(a, plot.it = F)$x)
}
write.csv(pheno, phenoOut, row.names = F)

if (all(pheno[, 1] == sqc.i[, 1]) == F){
  
  stop("ERROR: phenotype do not match sqc!")
}
sqc.i[, 12] <- ifelse(sqc.i[, 12] == "M", 1, 0)
write.csv(sqc.i[which(cnd), c(1:7)], sampleOut, row.names = F)
write.csv(sqc.i[, c(1,12, 27:36)], sqcOut, row.names = F)
