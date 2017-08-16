# DRD
script for the final exam of DNA RNA DYNAMICS course @ UNIBO

ASSIGNMENT:
1.
  Load raw data with minfi 
2.
  Perform the following quality checks: 
    QCplot 
    check the intensity of negative controls using shinyMethyl 
    check how many samples have >5 % of probes with a detection p-value greater than 0.01 
    check how many probes have > 1 % of samples with a detection p-value greater than 0.01 
    from this point on you can choose to continue using filtered or not filtered raw data 
3.
  Extract beta and M values and provide: 
    a density plot of their distributions, also dividing them according to type I and type II chemistry 
    a plot of mean methylation vs standard deviation 
    perform a PCA to check for batch effects 
4.
  Try different normalization methods and provide the density plots of mean values, standard deviation values and the boxplot of beta values for each sample. Choose a normalization method to perform differential analysis.
5.
  Perform differential analysis using a parametric and a not parametric test. If possible, apply the test on all the microarray probes. If not possible, select some chromosomes and apply the test on this smaller subset of probes. You can also choose to use covariates, for example BMI or Age. 
6.
  Apply multiple test correction and set a significant threshold of 0.05. How many probes do you identify as differentially methylated after Bonferroni correction? And how many after BH correction? 
7.
  Produce a Manhattan plot and a volcano plot of your data 
8.
  Produce an MDS plot on probes selected as differentially methylated. Produce an heatmap on probes selected as differentially methylated 
