# Genetic-Imputation
This code was developed as an in house solution to perform genetic imputation. Genetic samples can have many gaps in their data depending on the type of sequencing as well as the quality of sequencing. This method is used to predict a more complete human genome. 

This particular code was one of my first coding projects. I did not adequately use loops in my steps and never had time to go back and streamline my code.

This code does the following:
-Loads a reference genetic panel for each regional group.
-Initial data set up of:
    loading data
    removing SNPs with mendelian errors,
    ensure all strands are being read in the correct direction
    flipping SNPs if needed
    updating names of SNPs to match the reference panel's naming convention
    zipping the final imputation ready datasets
 - Phasing the reference panel and samples. This means we are estimating haplotypes in the dataset.
 - Converting the file to VCF format
 - Imputation of the sample using the reference panel using minimac3
 - Converting files to plink format for processing results
 - Creating text files to include covariates such as sex and phenotype
 - Adding covariate information into the plink file
 - Storing results
  
  These results were then used in a logistic regression using principle components to account for population genetics, sex, phenotype, Paternal ID, Maternal ID, Family ID, and all SNPs imputed with a high probablilty to ensure samples used are of best quality.
