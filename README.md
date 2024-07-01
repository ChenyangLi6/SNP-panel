Each Rmd file contains the code that was used to run the classification analysis from the .raw file for the corresponding data
There is a single summary file where we concatenated the results from each empirical classification called ClassificationResults_allComps.csv

We conducted analytical simulations where we fixed the allele frequency in pop1 = 1 and then generated GTs in pop2 to achieve a desired FST
The misclassification results from this are summarized in Simulate_FixAlleleFreqPop1_MissclassificationRate.csv

We also ran simulations in SLiM v4 that can be found in the TrousersModel_popSplit_SLiM_FSTOverTime_mainCoal2.slim file
The scripts classify_mergedVCF.R and makePlinkFiles.sh used to make .raw files and classify SLiM sims on cluster.
The summary file for SLiM classifications was called merged_classification_plotResults.txt and plotResults.R was used to plot the results. 
