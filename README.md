# DataAnalysis
This repository is for analyzing the [UKBiobank](http://www.ukbiobank.ac.uk/) data.

## Data Process
In this step, we select the samples and sort the order of phenotype data.

### STEP1: Combining the sqc data and fam data
The order of sqc data and fam data are the same. The sample size of the two data is 488377. <br>
The output of STEP 1: `QC sample: 488377`

### STEP2: Selecting samples from sqc data
Select the samples from combined sqc data.
The output of STEP 2:<br>
`Genotyping success: 487409 `<br/>
`White British ancestry subset: 408972` <br>
`Excess relatives: 188` <br>
`Sex chromosome aneuploidy: 652` <br>
`Used in PCA calculation: 407219` <br>
`Redacted: 14`<br>
`Samples Remaining: 377198`

### STEP3: Selecting and sorting samples of pehnotype data

### Reference
The introduction to datasets of UKBiobank: http://www.ukbiobank.ac.uk/wp-content/uploads/2017/07/ukb_genetic_file_description.txt

## Linear regression for each phenotypes
