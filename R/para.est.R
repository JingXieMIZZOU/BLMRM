
#' Hyperparameters Estimation
#'
#' A function to estimate hyperparameters, i.e., the scale and location parameters of prior distributions (Gamma) of the variances of random replicates effect and random SNP effect, as well as the parameters of the prior distribution (Gaussian) of fixed gene effect.This function tried to fit Generalized Linear Mixed Model to the data of each gene by \code{glmer} function in \code{lme4} and then filter out the genes with computational problems. The hyperparameter estimation is based on the genes without computational problems.
#'
#' @param data The raw data, need to be in the exact format with the sample data provided by this package, more specifically, a \code{data.frame} with its first 3 variables indicating gene ID, gene number and SNP number within a specific gene respectively. The following 2*Rep (Rep is the number of biological replicates) columns in the \code{data.frame} contain count data and the order must be consistent with the sample data. More details would be found in \code{help(mysample)}.
#' @param rep The number of biological replicates at each SNP.
#' @return A \code{list} contains the following:
#' \describe{
#' \item{para}{the estimation of hyperparameters.}
#' \item{index}{a vector contain index of genes without computational problems.}
#' \item{all}{a \code{data.frame} contains detailed intermediate results such as the p_values and estimated FDRs of likelihood ratio tests, the estimates of variance components etc.}
#' }
#' @import lme4
#' @import mvtnorm
#' @import MCMCpack
#' @import nloptr
#' @import dplyr
#' @import stats4
#' @import qvalue
#' @export




para.est<-function(data,rep){
    # data augumentation goes first
    # add two observations
    geneNum<-unique(data[,2])
    clean<- rep(NA,length(geneNum))
    P_SNP<- P_gene<- B4<- R4<- S4<- rep(NA,length(geneNum))
    for (i in geneNum){
        DATA<-GDD(i,data,rep=rep)
        D<-data.frame(SNP=DATA$SNP,NI=DATA$NI+2,YI=DATA$YI+1,Rep=DATA$Rep)
        tryCatch({
            print(i)
            mod2=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep),nAGQ = 0,
                       family = binomial, data = D)
            mod3=glmer(formula = cbind(YI, NI-YI)~ 0+(1 |Rep)+(1|SNP),
                       nAGQ = 0,family = binomial, data = D)
            mod4=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep)+(1|SNP),
                       nAGQ = 0, family = binomial, data =D)
            B4[i]<- as.vector(fixef(mod4))
            R4[i]<-lapply(VarCorr(mod4),'[[',1)$Rep
            S4[i]<-lapply(VarCorr(mod4),'[[',1)$SNP
            P_SNP[i]<-anova(mod2,mod4)$`Pr(>Chisq)`[2]
            P_gene[i]<- anova(mod3,mod4)$`Pr(>Chisq)`[2]
            clean[i]<-i
        }
        ,error=function(e){cat("error: ",conditionMessage(e), "\n")})
    }
    # organize resutls
    res<- data.frame(geneNum=na.omit(clean), P_gene=na.omit(P_gene),
                     P_SNP=na.omit(P_SNP), beta=na.omit(B4), Rep=na.omit(R4), SNP=na.omit(S4))
    res$q_gene<-qvalue(res$P_gene)$qvalues
    res$q_SNP<-qvalue(res$P_SNP)$qvalues
    set.seed(23457)
    tiebreaker <- sample(1:nrow(res), replace=TRUE)
    ## MOM to estimate paras of Inverse Gamma dist. for sigmaRsquare
    sigRf<-na.omit(res$Rep)
    uRs<-mean(sigRf)
    vRs<-var(sigRf)
    aRs<-uRs^2/vRs+2
    bRs<-uRs*(uRs^2/vRs+1)
    ## MOM to estimate paras of Inverse Gamma dist. for sigmaSsquareS
    res_SNP<-res[order(res$P_SNP,tiebreaker,decreasing = F),]
    sigSf<-na.omit(res_SNP[res_SNP$q_SNP<=0.05,]$SNP)
    uSs<-mean(sigSf)
    vSs<-var(sigSf)
    aSs<-uSs^2/vSs+2
    bSs<-uSs*(uSs^2/vSs+1)
    ## MLE to estimate paras of Gaussian prior for fixed gene effect
    res_gene<-res[order(res$P_gene,tiebreaker,decreasing = F),]
    betas<-res_gene[res_gene$q_gene<=0.05,]$beta
    betasf<-betas[!is.na(betas)]
    mlemus<- mean(betasf)
    mlesigmas<-sqrt(sum((betasf-mlemus)^2)/length(betasf))
    calibration<-data.frame(aRs=aRs,bRs=bRs,aSs=aSs, bSs=bSs,
                            mlemus=mlemus,mlesigmas=mlesigmas)
    return(list(para=calibration,index=na.omit(clean),all=res))
}
