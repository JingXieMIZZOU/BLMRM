#' Re-structurize the input data
#'
#' A function to convert the raw data into a structure which contains necessary information for analysis, i.e., SNPs, Replicates, counts  from maternal allel (YI) and total counts (NI). This structure is friendly to GLMM fitting in glmer() function in "lme4" package. This function can be easily modified to apply to other similar real data structures.
#'
#' @param data The raw data, need to be in the exact format with the sample data provided by this package, more specifically, a \code{data.frame} with its first 3 variables indicating gene ID, gene number and SNP number within a specific gene respectively. The following 2*Rep (Rep is the number of biological replicates) columns in the \code{data.frame} contain count data and the order must be consistent with the sample data. More details would be found in \code{help(mysample)}.
#' @param i The gene number of each gene. Take the sample dataset for instance, the gene number of the first gene in this dataset is 28.
#' @param rep The number of biological replicates at each SNP.
#' @return A \code{data.frame} in the structure that would be needed for analysis by BLRM method.
#' @export

GDD=function(i,data,rep){
    GeneDat=subset(data,data[,2]==i)
    S=length(unique(GeneDat[,3]))##number of SNP
    GD=matrix(0, ncol=4, nrow=rep*S)
    GD[,1]=kronecker(c(1:S),rep(1,rep))
    GD[,4]=kronecker(rep(1,S),c(1:rep))
    for(ss in 1:S){
        for(l in 1:rep){
            GD[(ss-1)*rep+l,2:3]=as.matrix(GeneDat[ss,((4+(l-1)*2):(5+(l-1)*2))])
        }
    }
    ##### Remove NA from Dataset
    GDNA= matrix(GD[complete.cases(GD),],ncol=4)
    GD=matrix(GDNA[(GDNA[,2]>0),],ncol=4)
    colnames(GD)=c("SNP","NI","YI","Rep")
    GD=data.frame(GD)
    GDV=data.frame(GD)
    GDV[,1]=as.factor(GDV[,1])
    GDV[,4]=as.factor(GDV[,4])
    return(GDV)
}


