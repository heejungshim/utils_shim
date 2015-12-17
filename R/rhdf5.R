
library(rhdf5)

##**************************************************************************
##  function to Get counts from hdf5
##  written by Ester Pantaleo and modified by Heejung Shim
##**************************************************************************
get.counts.h5 <- function(list_path, chr, locus.start, locus.end, list_name = NULL, print_message=FALSE){
    M <-  NULL
    for (h5file in list_path){
        if(print_message){
          print(paste0("Loading ", h5file))
          print(paste0("h5read(", h5file, ",", chr, ", index=list(", locus.start, ":", locus.end, ")"))
        }
        M        <- rbind(M, h5read(h5file, chr, index=list(locus.start:locus.end)))
    }
    row.names(M) = list_name
    return(M)
}


