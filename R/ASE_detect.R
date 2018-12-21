
#' ASE Detection
#'
#' A function to perform hypothesis tests.
#'
#' @param data The raw data, need to be in the exact format with the sample data provided by this package, more specifically, a \code{data.frame} with its first 3 variables indicating gene ID, gene number and SNP number within a specific gene respectively. The following 2*Rep (Rep is the number of biological replicates) columns in the \code{data.frame} contain count data and the order must be consistent with the sample data. More details would be found in \code{help(mysample)}.
#' @param clean_index a vector contain index of genes without computational problems which would be returned by \code{par.est} function.
#' @param paras hyperparameters estimates returned by \code{par.est} function.
#' @param rep The number of biological replicates at each SNP.
#' @param fdr The desired False Discovery Rate (FDR), default value was set as 0.05.
#'
#' @return A \code{list} contains the following:
#' \describe{
#' \item{geneEffect}{A \code{data.frame} contains gene IDs and gene numbers of the genes with significant ASE and their corresponding posterior probabilities and estimated FDRs.}
#' \item{SNPEffect}{A \code{data.frame} contains gene IDs and gene numbers of the genes with significant ASE variation across SNPs and their corresponding posterior probabilities and estimated FDRs.}
#' \item{GSEffect}{A \code{data.frame} contains gene IDs and gene numbers of the genes exhibiting both ASE gene effect and ASE variation across SNPs and their corresponding posterior probabilities and estimated FDRs.}
#' \item{logPPs}{a \code{data.frame} contains detailed intermediate results, i.e., the posterior probabilities of the 4 models for all the genes.}
#' }
#' @import lme4
#' @import mvtnorm
#' @import MCMCpack
#' @import Matrix
#' @import nloptr
#' @import dplyr
#' @import stats4
#' @import qvalue
#' @import IDPmisc
#' @export


detection<-function(data,clean_index,paras,rep,fdr=0.05){

    ## realPPS is the log of PP before standardize
    temp<-lapply(clean_index,PPsFUN,par=paras,data=data,rep=rep)
    realPP<-do.call(rbind, lapply(temp, head, 1))

    ## normalized PPs
    realPPS<-exp(realPP[,2:5])
    realPPSS<-matrix(nrow=nrow(realPP),ncol=6)
    realPPSS<-as.data.frame(realPPSS)
    realPPSS[,1]<-realPP[,1]
    realPPSS[,2]<-unlist(lapply(as.vector(realPP[,1]),function(x) unique(data[data[,2]==x,1])))
    realPPSS[,3]<-realPPS[,1]/(apply(realPPS,1,sum))
    realPPSS[,4]<-realPPS[,2]/(apply(realPPS,1,sum))
    realPPSS[,5]<-realPPS[,3]/(apply(realPPS,1,sum))
    realPPSS[,6]<-realPPS[,4]/(apply(realPPS,1,sum))
    colnames(realPPSS)<-c("GeneNum","GeneID","PP1","PP2","PP3","PP4")

    # re-identification
    # gene effect
    real_gene<-data.frame(geneNum=realPPSS$GeneNum,
                          geneID=realPPSS$GeneID,
                          PP=realPPSS$PP2+realPPSS$PP4,
                          FPP=1-(realPPSS$PP2+realPPSS$PP4))
    real_gene<-real_gene[order(real_gene$PP,decreasing = TRUE),]
    real_gene$FDR<-cummean(real_gene$FPP)
    real_gene_res<-real_gene[real_gene$FDR<=fdr,]

    # SNP effect
    # re-identification
    real_SNP<-data.frame(geneNum=realPPSS$GeneNum,
                         geneID=realPPSS$GeneID,
                         PP=realPPSS$PP3+realPPSS$PP4,
                         FPP=1-(realPPSS$PP3+realPPSS$PP4))
    real_SNP<-real_SNP[order(real_SNP$PP,decreasing = TRUE),]
    real_SNP$FDR<-cummean(real_SNP$FPP)
    real_SNP_res<-real_SNP[real_SNP$FDR<=fdr,]

    # Gene & SNP effect
    real_GS<-data.frame(geneNum=realPPSS$GeneNum,
                        geneID=realPPSS$GeneID,
                        PP=realPPSS$PP4,
                        FPP=1-(realPPSS$PP4))
    real_GS<-real_GS[order(real_GS$PP,decreasing=TRUE),]
    real_GS$FDR<-cummean(real_GS$FPP)
    real_GS_res<-real_GS[real_GS$FDR<=fdr,]

    return(list(GeneEffect=na.omit(real_gene_res),
                SNPEffect=na.omit(real_SNP_res),
                GSEffect=na.omit(real_GS_res),
                logPPs=realPPSS
    ))
}
