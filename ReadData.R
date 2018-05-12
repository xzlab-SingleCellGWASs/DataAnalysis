###This code is for transfering rds data to the normal counts data
library(SingleCellExperiment)
library(dplyr)

###Read the RDS data
data=readRDS("~/darmanis.rds")

### just type data, you will see there are information in the assay like : assays(3): counts normcounts logcounts
data

###if the assay contains counts, then obtain the count data from RDSdata
raw_counts=as.matrix(counts(data))
###if the assay contains normcounts, then obtain the count data from RDSdata
norm_counts=as.matrix(normcounts(data))
###if the assay contains logcounts, then obtain the count data from RDSdata
log_counts=as.matrix(logcounts(data))

### obtain the gene name
genename=rownames(data)

### cell type information
colData=colData(data)
### you can list colnames(colData) there are some infromation in it.
### depends on the rds format, some have the cell_type infromation stored in the cell_type1 or cell_type2 column
### for example
celltype=colData$cell_type1
