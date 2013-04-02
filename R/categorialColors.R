.categorialColors <- function(n){

  if(n<=17){
    categorial_17 <- c(30,36,51,5,10,24,33,7,48,96,117,122,26,142,190,12,1)
    col.n <- categorial_17[1:n]
    return(colors()[col.n])
  }else{
    stop("Number of colors must be smaller or equal than 17.")
  }
  
}

