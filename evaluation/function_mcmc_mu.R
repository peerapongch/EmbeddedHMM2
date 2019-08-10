mci_mu <- function(sample,T,dim){
  mcmc_mu <- matrix(0,nrow=T,ncol=dim)
  for(t in 1:T){
    mcmc_mu[t,] <- apply(sample[,t,],MARGIN=2,FUN=mean)
  }
  return(mcmc_mu)
}
