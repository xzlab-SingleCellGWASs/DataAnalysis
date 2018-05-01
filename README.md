# DataAnalysis
This repository is for analyzing the [UKBiobank](http://www.ukbiobank.ac.uk/) data.

## Preparing datasets
In this section, we select the samples and sort the phenotype data.

### STEP1: Combining the sqc data and fam data
The order of sqc data and fam data are the same. The sample size of the two data is 488377. <br>
The output of STEP 1: `QC sample: 488377`

### STEP2: Selecting samples from new sqc data
Select the samples from combined sqc data.<br>
The output of STEP 2:<br>
`Genotyping success: 487409 `<br/>
`White British ancestry subset: 408972` <br>
`Excess relatives: 188` <br>
`Sex chromosome aneuploidy: 652` <br>
`Used in PCA calculation: 407219` <br>
`Redacted: 14`<br>
`Samples Remaining: 377198`

### STEP3: Selecting and sorting samples of phenotype data
Sort phenotype data as the order of sqc data because order of the genotype data is same to the sqc data. 

### Reference
The introduction to datasets of UKBiobank: http://www.ukbiobank.ac.uk/wp-content/uploads/2017/07/ukb_genetic_file_description.txt


## Analysis UKBiobank data by cross validation
In this section, we get the summary data.

### STEP1: Setting parameters
`minMAF = 1e-3 `<br/>
`minINFO = 0.8`<br/>
`minHW = 1e-10`<br/>
`callingRate = 0.95`<br/>
`prop = 0.8;`<br/>

### STEP2: Loading and processing data
Load the data from the first section, including phenotype, sqc and sqcNA.
Get the index of each samples. The selection standards are sqc, phenotype and cross validation.

### STEP3: Getting summary data
Output the summary data of all selected SNPs.
Get the infomation of SNPs < 1e-3 and < 1e-8 for the further analysis.

