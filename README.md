# CNPFP
Protein Function Prediction
##Introduction

CNPFP is a network based method for protein function prediction. The codes and data are available for reproducibility.

##Datasets

CNPFP was tested using 2 datasets

*Yeast Cell Cycle;
*Yeast Protein Protein Interaction
The annotation dataset can be downloaded from http://frasca.di.unimi.it/cosnetdata/GO.ann.yeast.28.03.13.3_300.txt

To load the datasets please type:
source("(file location)/YPPINetworkFunctions.txt");
source("(file location)/YPPINetworkFunctions1.txt");
read.csv(file = "(file location)/YeastProtein_NodeNames.csv", header=F)
read.csv("(file location)/YeastProtein_ADJ1.csv", header=F)
read.csv(file ="(file location)/YEASTCellCycle4000.csv",header=T, row.names=1)

To run CNPFP execute 'CNPFP.r'

#Evaluation 

auroc
accuracy
F-Score
