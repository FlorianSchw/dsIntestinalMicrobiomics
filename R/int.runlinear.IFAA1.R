#### Needs proper Documentation as an internal function






int.runlinear.IFAA1=function(
    x,
    y,
    nPredics,
    fwerRate=0.25,
    adjust_method="fdr"
){

  results=list()


  bootResu <- int.lm_sparse(x = x, y = y)

  bootResu[[9]] <- nPredics
  bootResu[[10]] <- fwerRate






  return(bootResu)
}
