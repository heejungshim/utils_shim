
##**************************************************************************
##  Function to make given colors transparent. 
##  This function is copied from http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color 
##  Usage : makeTransparent(c("red","blue"), 200) or makeTransparent("red", 200)
##'
##' @param sameColor; a vector of colors
##' @param alpha; scaler, value between 0 and 255
##' @return a vector of colors
##**************************************************************************
makeTransparent <- function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}




