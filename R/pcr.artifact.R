## This file contains useful functions for analysis.
## Sometime I copied original functions and modified them. Then, I specified original sources. 
##
## Copyright (C) 2014 Heejung Shim
##



##' 'get.pcr.artifacts.posi' checks if any sample has pcr artifacts in a given position and returns index of samples (a vector) that are identified to have pcr artifacts in the given position. 
##'
##' Here, we consider a window of left and right `win.half.size' bps around a given position. If read count at the given positon is more than `prop.thresh' of reads in the window, we consider read count at the given position is due to pcr artifacts.
##' 
##' 
##' @param posi a candidate position for pcr artifacts in data 
##' @param data a matrix of num of samples by site size; original count data of all samples in a given site
##' @param win.half.size default=50; a half of window size
##' @param prop.thresh default=0.9;  
##' @return pcr.sample a vector of samples that are identified to have pcr artifacts in a given position. 
get.pcr.artifacts.posi <- function(posi, data, win.half.size = 50, prop.thresh = 0.9){
   
    st.win = max(1, posi - win.half.size)        
    en.win = st.win + win.half.size*2
    if(en.win > dim(data)[2]){
        en.win = dim(data)[2]
        st.win = en.win - win.half.size*2
    }
    prop = data[,posi]/apply(data[,st.win:en.win], 1, sum)
    pcr.sample = which(prop > prop.thresh)
    return(pcr.sample)
}



##' 'remove.pcr.artifacts' identifies positions in data where at least one sample has pcr artifacts (using the function `get.pcr.artifacts.posi') and replace read counts of all samples at that position with max(1, average of read counts of sample without pcr artifacts).
##'
##'
##' 
##' 
##' @param data a matrix of num of samples by site size; original count data of all samples in a given site
##' @param win.half.size default=50; argument to get.pcr.artifacts.posi
##' @param prop.thresh default=0.9; argument to get.pcr.artifacts.posi
##' @return a list of data and posi.with.pcr.artifacts; data contains pcr artifacts-removed data; posi.with.pcr.artifacts contains a list of positions where at least one sample has pcr.artifacts. 
remove.pcr.artifacts <- function(data, win.half.size = 50, prop.thresh = 0.9){

    # only consider positions with at least 2 reads as a candidate position 
    max.val = apply(data, 2, max)
    candidate.posi = which(max.val > 1)

    num.sam = dim(data)[1]
    len = length(candidate.posi)
    if(len > 0){
        # for each candidate position, get which samples have pcr artifacts
        pcr.posi.list = lapply(candidate.posi, get.pcr.artifacts.posi, data = data, win.half.size = win.half.size, prop.thresh = prop.thresh)
        pcr.posi = which(sapply(pcr.posi.list, length) > 0)
        len.pcr.posi = length(pcr.posi)
        if(len.pcr.posi > 0){
            # for positions where at least one sample has pcr artifacts, replace read counts of all samples with max(1, average of read counts of sample without pcr artifacts)
            for(p in 1:len.pcr.posi){
                pcr.sam = pcr.posi.list[[pcr.posi[p]]]
                ix = candidate.posi[pcr.posi[p]]
                if(length(pcr.sam) == num.sam){
                    data[,ix] = 1
                }else{
                    data[, ix] = max(1, ceiling(mean(data[-pcr.sam,ix])))
                }
            }
        }
    }else{
        len.pcr.posi=0
    }
    
    if(len.pcr.posi > 0){
        return(list(data=data, posi.with.pcr.artifacts = candidate.posi[pcr.posi]))
    }else{
        return(list(data=data, posi.with.pcr.artifacts = NULL))
    }
}


##' 'remove.pcr.artifacts.in.known.posi' uses pre-identified positions where at least one sample has pcr artifacts and replace read counts of all samples at that position with max(1, average of read counts of sample without pcr artifacts).
##'
##'
##' 
##' 
##' @param data a matrix of num of samples by site size; original count data of all samples in a given site
##' @param known.pcr.posi a vector of positions 
##' @param win.half.size default=50; argument to get.pcr.artifacts.posi
##' @param prop.thresh default=0.9; argument to get.pcr.artifacts.posi
##' @return a list of data and posi.with.pcr.artifacts; data contains pcr artifacts-removed data; posi.with.pcr.artifacts contains a list of positions where at least one sample has pcr.artifacts. 
remove.pcr.artifacts.in.known.posi <- function(data, known.pcr.posi, win.half.size = 50, prop.thresh = 0.9){

  if(length(known.pcr.posi) > 0){
    num.sam = dim(data)[2]
  
    pcr.posi.list = lapply(known.pcr.posi, get.pcr.artifacts.posi, data = data, win.half.size = win.half.size, prop.thresh = prop.thresh)

    len.pcr.posi = length(pcr.posi.list)
    for(p in 1:len.pcr.posi){
      pcr.sam = pcr.posi.list[[p]]
      ix = known.pcr.posi[p]
      if(length(pcr.sam) == num.sam){
        data[,ix] = 1
      }else{
        data[, ix] = max(1, ceiling(mean(data[-pcr.sam,ix])))
      }
    }
  }
  return(list(data=data))
}


