###This program aims to processing dataset, including phenotye and QC. 

library(data.table)
library(dplyr)
library(pipeR)
library(optparse)

#####Input and output files 
###Input
##sqc data: ukb_sqc_v2.txt
##sqc head data: ukb_sqc_v2_head.txt (in the software)
##fam data: ukb30186_cal_chrXX_v2_s488363.fam (Any fam file is OK)
##phenotype: ukb10683.csv (RData)
##phenotype list: lable of analyzed phenotypes
###Output
##output path
args_list = list(
            make_option(c("-s", "--sqc"), type="character", default=NULL, 
                        help="INPUT: sqc data", metavar="character"),
            make_option(c("-c", "--colSqc"), type="character", default=NULL, 
                        help="INPUT: col names of sqc data", metavar="character"),
            make_option(c("-f", "--fam"), type="character", default=NULL, 
                        help="INPUT: fam file", metavar="character"),
            make_option(c("-p", "--pheno"), type="character", default=NULL, 
                        help="INPUT: phenotype file", metavar="character"),
            make_option(c("-l", "--phenoList"), type="character", default=NULL, 
                        help="INPUT: phenotype for analyzing", metavar="character"),
            make_option(c("-o", "--output"), type="character", default=NULL, 
                        help="Output path", metavar="character")
) 

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

#########STEP1: Combine the sqc and fam
print("STEP1: Combine the sqc and fam...")
#File input
sqc <- fread(opt$sqc, stringsAsFactors = F)
sqc <- data.frame(sqc)
header <- read.table(opt$colSqc, stringsAsFactors = F)
header <- data.frame(header)
fam <- fread(opt$fam, stringsAsFactors = F)
fam <- data.frame(fam)

#
if (dim(fam)[1] != dim(sqc)[1]){
 
  stop("ERROR: The number of rows is not the same!")
}else{
  
  sqc <- cbind.data.frame(fam[1], sqc)
  names(sqc) <- c("eid", header[, 1])
  cat("QC sample:", dim(sqc)[1], "\n")
}

#########STEP2: Select the samples

##in.Phasing.Input.chr1_22==1
##in.white.British.ancestry.subset==1
##used.in.pca.calculation==1
##excess.relatives==0
##putative.sex.chromosome.aneuploidy==0
##eid > 0
print("STEP2: Select the samples...")
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
write.csv(as.numeric(cnd), paste0(opt$output, "cndOut.csv"), row.names = F)
write.csv(sqc.i[which(cnd), c(1:7)], paste0(opt$output, "sampleInfo.csv"), row.names = F)

#########STEP3: phenotype data
print("STEP3: Process phenotype data...")
load(opt$pheno)
label <- read.table(opt$phenoList, header = F)

#Select and sort the phenotype data
pheno <- ukb[which(ukb$eid %in% eid.i), ]
idx <- match(eid.i, pheno[, 1])
pheno <- pheno[idx, ]

if (all(label[, 1] %in% colnames(pheno)) == F){
  
  stop("ERROR: The labels of phenotype are wrong!")
}
pheno <- pheno[, c(1, which(colnames(pheno) %in% label[, 1]))]
pheno[which(is.na(pheno[, 1])), 1] <- c(-1:-14)


#########STEP4: regress between pheno and sqc
print("STEP4: Regress between pheno and sqc...")
#Check the size of sqc and pheno
if (all(pheno[, 1] == sqc.i[, 1]) == F){
  
  stop("ERROR: phenotype do not match sqc!")
}
#Regress
sqc.i[, 12] <- ifelse(sqc.i[, 12] == "M", 1, 0)
if (nrow(label) == 1){
 
  pheno[, 2] <- qqnorm(pheno[, 2], plot.it = F)$x
  phenoSqc <- data.frame(cbind(pheno[, 2], sqc.i[, c(12, 27:36)]))
  colnames(phenoSqc)[1] <- "y"
  resid <- resid(lm(y ~ ., data = phenoSqc)) %>>%
            {qqnorm(., plot.it = F)$x}
  write.csv(resid, file = paste0(opt$output, "phenoResid.csv"),
            row.names = F)
}else{
  
  pheno[, c(2: ncol(pheno))] <- apply(pheno[, c(2: ncol(pheno))], 2, function(a) 
                           qqnorm(a, plot.it = F)$x)
  for (i in 2: c(ncol(pheno))){
    
    phenoSqc <- data.frame(cbind(pheno[, i], sqc.i[, c(12, 27:36)]))
    colnames(phenoSqc)[1] <- "y"   
    resid <- resid(lm(y ~ ., data = phenoSqc)) %>>%
    {qqnorm(., plot.it = F)$x}  
    if (i < 10){
      
      pheno.code <- paste0("0", i)
      write.csv(resid, file = paste0(opt$output, "pheno", pheno.code, "Resid.csv"),
                row.names = F)
    }else{
      
      write.csv(resid, file = paste0(opt$output, "pheno", i, "Resid.csv"), 
                row.names = F)
    }
  }
}
print("ALL DONE!")