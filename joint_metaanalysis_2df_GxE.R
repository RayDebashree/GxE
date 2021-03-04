###------------- Code for implementing joint meta-analysis of 2-df GxE effect -----------------

## Code developed by: Debashree Ray, PhD, MStat
## Contact info		: dray@jhu.edu
## Last modified		: Aug-10-2020
## Based on Manning et al (2011) Genetic Epidemiology paper (PMID: 21181894)
## Performs meta-analysis and then performs the 2-df chisq test of H0: ß_G=0, ß_GE=0

## Cite the following articles if this function is used:
##	Manning, A.K., LaValley, M., Liu, C.T., ... , Dupuis, J. (2011) "Meta-analysis of gene-environment interaction: 
##	joint estimation of SNP and SNP x environment regression coefficients". Genet Epidemiol 35(1), 11-18, https://doi.org/10.1002/gepi.20546
## Zhang, W., Venkataraghavan, S., Hetmanski, J.B., ... , Ray, D., Beaty, T.H. (2021+, under review) "Detecting 
##	gene-environment interaction for maternal exposures using case-parent trios ascertained through a case with nonsyndromic orofacial cleft".

message("==========================================")
message("     Joint 2-df GxE meta-analysis")
message("==========================================")
message("")

JMA2df <- function(SNP.coef, GXE.coef, SNP.se, GXE.se, SNP.GXE.cov=NULL, trio.package.input=FALSE){
    # checks
    k <- length(SNP.coef)
    if(length(GXE.coef)!=k | length(SNP.se)!=k | length(GXE.se)!=k)
	stop("Inputs SNP.coef, GXE.coef, SNP.se, GXE.se should be vectors of same length k, where k=number of studies to meta-analyze.")	
    if(!(trio.package.input) & is.null(SNP.GXE.cov))
	stop("Input SNP.GXE.cov (vector of covariance between SNP.coef and GXE.coef) should be provided if the inputs are NOT from R trio package.")
    if(!is.null(SNP.GXE.cov) & length(SNP.GXE.cov)!=k)
	stop("Inputs SNP.coef, GXE.coef, SNP.se, GXE.se, and SNP.GXE.cov should be vectors of same length k, where k=number of studies to meta-analyze.")
	if(k<2) stop("You need to provide information on at least 2 studies - else there's nothing to meta-analyze here.")
    if(trio.package.input)
	SNP.GXE.cov <- (-1)*(SNP.se^2)

    # define vectors/matrices
    b <- matrix(unlist(lapply(1:k, function(i) c(SNP.coef[i],GXE.coef[i]))), ncol=1)
    Covb <- matrix(0, nrow=2*k, ncol=2*k)
    diag(Covb) <- unlist(lapply(1:k, function(i) c(SNP.se[i]^2,GXE.se[i]^2) ))
    for(i in 1:k){
	Covb[(2*i-1),(2*i)] <- Covb[(2*i),(2*i-1)] <- SNP.GXE.cov[i]
    }
    W <- matrix(0, nrow=2*k, ncol=2)
    W[seq(1,2*k,2),1] <- W[seq(2,2*k,2),2] <- 1

    # meta-analyzed estimates from generalized least squares
    mat.1 <- t(W) %*% solve(Covb)
    Cov.beta.hat.inv <- mat.1 %*% W
    Cov.beta.hat <- solve(Cov.beta.hat.inv)
    beta.hat <- Cov.beta.hat %*% mat.1 %*% b

    # 2-df Joint meta-analysis test of H0: ß_G=0, ß_GE=0
    stat <- c(t(beta.hat) %*% Cov.beta.hat.inv %*% beta.hat)
    pval <- pchisq(stat, df=2, ncp=0, lower.tail=F)

    return(list(SNP.coef.JMA=beta.hat[1], GXE.coef.JMA=beta.hat[2], SNP.se.JMA=sqrt(Cov.beta.hat[1,1]), GXE.se.JMA=sqrt(Cov.beta.hat[2,2]), 
SNP.GXE.cov.JMA=Cov.beta.hat[1,2], wald2df.stat.JMA=stat, wald2df.pval.JMA=pval))
}
