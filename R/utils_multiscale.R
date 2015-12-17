## This file contains useful functions for multiscale analysis 
## Sometime I copied original functions and modified them. Then, I specified original sources. 
##
## Copyright (C) 2014 Heejung Shim
##


##' 'compute.proportions.given.effect.size' takes effect sizes over a retion (as a ratio)
##' and sum of proportions (in binomial sampling) for two groups as input. And it computes
##' proportion parameters in binomial sampling for two groups for the given effect sizes
##' (as a ratio). Specifically, let $p_1 = p$ and $p_2 = p \times effect.ratio$
##' where $\log(effect.ratio)$ is effec.size in log space and $effect.ratio = 1$ means
##' there is no effect. Use $p$ which satisfies $p_1 + p_2 = sum.prop$.
##' So $p_1 = p = sum.prop\frac{1}{1+effect.ratio}$ and $p_2 = sum.propr\frac{effect.ratio}{1+effect.ratio}$.
##'
##'
##' for example,
##' sum.prop = rep(2/70, 1024)
##' effect.ratio = rep(1, 1024)
##' effect.ratio[100:200] = 1/2
##' res = compute.proportions.given.effect.size(sum.prop = sum.prop, effect.ratio = effect.ratio)
##' res$prop0
##' res$prop1
##' 
##' @param sum.prop a vector of sum of proportions (in binomial sampling) for two groups over a region
##' @param effect.ratio a vector of effect sizes over a region as a ratio; default value = NULL; if effect.ratio == NULL, we assign the same proportion to two groups (no effect).
##' @return a list of prop0 and prop1; prop0 (prop1) is a vector of proportions (in binomial sampling) for group0 (group1) over a region 
compute.proportions.given.effect.size <- function(sum.prop, effect.ratio=NULL) {

  if(is.null(effect.ratio)){
    effect.ratio = rep(1, length(sum.prop))
  }
  prop0 = sum.prop/(effect.ratio + 1)
  prop1 = sum.prop/(effect.ratio + 1)*effect.ratio
  
  return(params = list(prop0 = prop0, prop1 = prop1))
}



##' 'estBetaParams' takes mean and variance for beta distribution and returns
##' alpha and beta parameters in beta distribution.
##'
##'
##' I copied the original version from http://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
##' Then I modified to handle cases when somme elements of mu are 0 or 1 (then it returns NA if mu <= 0 or >=1) or when computed alpha and beta <= 0 (then it returns NA)
##'
##' @param mu.orig a vector of mean for beta distribution
##' @param var a vector (or scalar) of variance for beta distribution
##' @return a list of alpha and beta 
estBetaParams <- function(mu.orig, var) {

    del.ix = ((mu.orig <= 0) | (mu.orig >= 1))
    if(sum(del.ix) > 0){
        mu = mu.orig[!del.ix]
    }else{
        mu = mu.orig
    }
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)

    if(sum(del.ix) > 0){
        alpha.orig = beta.orig = rep(NA, length(mu.orig))
        alpha.orig[!del.ix] = alpha
        beta.orig[!del.ix] = beta
    }else{
        alpha.orig = alpha
        beta.orig = beta
    }

    ## handle non-positive alpha or beta 
    invalid.para = which((alpha.orig <= 0) | (beta.orig <= 0))
    if(length(invalid.para) > 0){
        alpha.orig[invalid.para] = rep(NA, length(invalid.para))
        beta.orig[invalid.para] = rep(NA, length(invalid.para))
    }
    
    return(params = list(alpha = alpha.orig, beta = beta.orig))
}





##' 'sample.from.Binomial.with.Overdispersion' simulates binomial samples with/without
##' over dispersion. 
##'
##' For a given overdispersion parameter, computed parameters for beta distribution can be invalid (e.g., mu.sig are 0 or 1). Then we sample read without overdispersion for those positions.
##' 
##' @param num.sam number of samples to be sampled
##' @param total.count a vector of non-negative counts;
##' @param mu.sig a vector of probabilities (we allow 0 or 1 as probablity)
##' @param over.dispersion if over.dispersion == NULL, simulate data from binomial. If over.dispersion is provided, simulate data binomial with over.dispersion.
##' @return a matrix of num.sam by L (length of total.count) containing simulated data. 
sample.from.Binomial.with.Overdispersion <- function(num.sam, total.count, mu.sig, over.dispersion=NULL){

    invalid.entry = ((mu.sig < 0) | (mu.sig > 1))
    if(sum(invalid.entry) > 0){ stop("ERROR, mu.sig have some values outside of valid range [0, 1]")}
                                
     if(is.null(over.dispersion)){
        return(matrix(data=rbinom(length(mu.sig)*num.sam, total.count, mu.sig), nr = num.sam, byrow = TRUE))
    }else{

        
        final.dat = matrix(data=NA, nr = num.sam, nc = length(mu.sig))
        
        # get alpha and beta
        resBeta = estBetaParams(mu.sig, over.dispersion)
        alpha = resBeta$alpha
        beta = resBeta$beta

        # for valid alpha and beta, sample data 
        del.ix = is.na(alpha)
        p.sig = rbeta(sum(!del.ix)*num.sam, alpha[!del.ix], beta[!del.ix])
        dat.V = rbinom(sum(!del.ix)*num.sam, total.count[!del.ix], p.sig) 
        final.dat[,!del.ix] = matrix(data=dat.V, nr = num.sam, byrow = TRUE)

        # for invalid alpha and beta, sample without over dispersion
        if(sum(del.ix) > 0){
            dat.IV = matrix(data=rbinom(sum(del.ix)*num.sam, total.count[del.ix], mu.sig[del.ix]), nr = num.sam, byrow = TRUE)
            final.dat[,del.ix] = matrix(data=dat.IV, nr = num.sam, byrow = TRUE)
        }

        return(final.dat)

    }

}
 



