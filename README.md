# DataAnalysis
This repository is for analyzing the [UKBiobank](http://www.ukbiobank.ac.uk/) data.

## Section 1: Preparing datasets
In this section, we select the samples and sort the phenotype data.

### STEP 1: Combining the sqc data and fam data
The order of sqc data and fam data are the same. The sample size of the two data is 488377. <br>
The output of STEP 1: `QC sample: 488377`

### STEP 2: Selecting samples from new sqc data
Select the samples from combined sqc data.<br>
The output of STEP 2:<br>
`Genotyping success: 487409 `<br/>
`White British ancestry subset: 408972` <br>
`Excess relatives: 188` <br>
`Sex chromosome aneuploidy: 652` <br>
`Used in PCA calculation: 407219` <br>
`Redacted: 14`<br>
`Samples Remaining: 377198`

### STEP 3: Selecting and sorting samples of phenotype data
Sort phenotype data as the order of sqc data because order of the genotype data is same to the sqc data. 

### Reference
The introduction to datasets of UKBiobank: http://www.ukbiobank.ac.uk/wp-content/uploads/2017/07/ukb_genetic_file_description.txt


## Section 2: Analyzing UKBiobank data by cross validation
In this section, we get the summary data.

### STEP 1: Setting parameters
`minMAF = 1e-3 `<br/>
`minINFO = 0.8`<br/>
`minHW = 1e-10`<br/>
`callingRate = 0.95`<br/>
`prop = 0.8;`<br/>

### STEP 2: Loading and processing data
Load the data from the first section, including phenotype, sqc and sqcNA.
Get the index of each samples. The selection standards are sqc, phenotype and cross validation.

### STEP 3: Getting summary data
Using the function `summ`, we can easily get the summary data of all selected SNPs from the bgen format. <br>
```
summ -maf mad_num -info info_num -hwe hwe_num -call call_num \
     -thread thread_num -prop prop_num -seed sedd_num \
     -pheno pheno_file -sqc sqc_file -bgen bgen_file \
     -outpath out_path -outfile out_file -chr chr_num -cv 0
```
- -maf, -info, -hwe and -call: the minimum of MAF, information, hwe and calling rate.
- -thread: the thread to parallel.
- -prop and -seed: the proportion of training data and seed.
- -pheno: phenotype data (csv format).
- -sqc: sqc index data (csv format).
- -bgen: bgen data.
- -outpath and -output: outpath and outfile (txt format)
- -chr: chromosome number.
- -cv: cv number.

