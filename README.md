# BLRM: A Bayesian Logistic Regression Model for Genome-Wide Detection of Allele-Specific Gene Expression

This R package provide a powerful and flexiable method (BLRM) to detect ASE genes and detect ASE variation within genes simultaneously while maintaining low computational requirements.

## Getting Started
There are basically three steps to conduct ASE detection with BLRM, i.e., read in raw data, estimate hyperparameters, and test hypotheses.Typical ASE analysis with BLRM can be performed by running below code, see introduction.pdf inside vignettes folder for details.

```
library("BLRM")
rawdata<-read.csv(file="YourRawdata.csv")
hyperparas<- para.est(data=rawdata,rep=R)
res<- detection(data=rawdata,clean_index=hyperparas$index,paras=hyperparas$para,rep=R,fdr=0.05)
list.ASEgene<-res$GeneEffect
list.SNPvariation<-res$SNPEffect
list.ASE.SNP<- res$GSEffect
```

## Prerequisites
R (>= 2.10)

## Installing

```
install.packages("devtools")
devtools::install_github("JingXieMIZZOU/BLRM")
```


## Authors
 Jing Xie (jx5fd@mail.missouri.edu)

## License
This project is licensed under the GPL License.

