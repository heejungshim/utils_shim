## This file contains useful functions for analysis.
## Sometime I copied original functions and modified them. Then, I specified original sources. 
##
## Copyright (C) 2014 Heejung Shim
##

##' 'get.pval.from.empirical.null.dist.discrete' takes empirical distribution of statistic
##' under the null (statistic.null) and a series of observed test statistics (statistic.alt),
##' and returns a series of p-values corresponding to the observed test statistics.
##' When test statistic is not completely continuous, there are problems for the qvalue
##' package which uses the insignificant p values to estimate the proportion of nulls.
##' Here, this function uses randomization to produce continuous p-values that are uniformly
##' distributed under the null.
##'
##' Let t be observed test statistic. p-value is defined by P(T >= t | H0), but this is not
##' uniformly distributed unless test statistic is continuous. We randomize them using
##' P(T > t | H0) + U*P(T = t | H0), where U ~ uniform(0,1). The proposed p-value is
##' between P(T > t | H0) and P(T > t | H0) + P(T = t | H0). This p-value can be empirically
##' computed by #[T > t] + U*(#[T = t] + 1) / (total number of statistic under the null + 1).
##' Here "+1" in the denominator is due to observed t. 
##'
##' @param statistic.null a vector of test statistics under the null
##' @param statistic.alt a vector of observed test statistics
##' @param big.sig bool indicating whether bigger statistic is more significant.
##' @return a vector of p-values corresponding to the observed test statistics.
get.pval.from.empirical.null.dist.discrete <- function(statistic.null, statistic.alt, big.sig = TRUE){
    
    numNulltests = length(statistic.null)
    if(big.sig){
        numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null > x))}, statistic.null = statistic.null)
      }else{
        numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null < x))}, statistic.null = statistic.null)
    }
    numEqual = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null == x))}, statistic.null = statistic.null)
    Uval = runif(length(statistic.alt))
    
    pval.list = (numSig + Uval*(numEqual + 1))/(numNulltests + 1)
    return(pval.list)
}

