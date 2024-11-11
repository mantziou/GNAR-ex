################################################################
######## CORELLATIONS BETWEEN GDP AND PAYMENTS TS ########
################################################################

pay_gdp_corr_crossind_egerem <- function(prep_dat,no_nodes,no_edges,thres,growth_rates=TRUE){
  # prep_dat: data verion i, data frame with edge and nodal time series
  # no_nodes: number of nodes in version i
  # no_edges: number of edges in version i
  # thres: threshold for Pearson's correlation 
  # return: payments (edges) to be removed according to Pearson's threshold
  
  if (growth_rates){
    data <- growthrates(prep_dat)
  }
  cormat <-  matrix(,ncol=no_nodes,nrow = no_edges)
  for (i in 1:no_edges){
    ind_col <- 1
    for (j in (no_edges+1):(no_nodes+no_edges)){
      cormat[i,ind_col] <- cor(data[i,],data[j,],method = "pearson")
      ind_col <- ind_col+1
    }
  }
  outlierspos <- which((cormat<(-thres))|(cormat>thres),arr.ind = TRUE)
  paymrem <- unique(outlierspos[,1])
  return(list(paymrem=paymrem,cormat=cormat))
}

