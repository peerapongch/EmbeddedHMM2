ACTime <- function(mcmcs,T,dim,lag.max,tps=1){
  # print(tps)
  # calculate overall mean 
  overall_mean <- matrix(0,nrow=T,ncol=dim)
  N <- dim(mcmcs[[1]]$X_sample)[1]
  remove <- round(N*0.1)
  N <- N-remove # reassign N  
  divisor <- N*length(mcmcs)
  # ac_time <- array(0,dim=c(length(mcmcs),T,dim))
  ac_time <- matrix(0,nrow=T,ncol=dim)
  pb <- txtProgressBar(min=0,max=length(mcmcs)*T*dim,style=3); prog <- 0
  for(t in 1:T){
    for(j in 1:dim){
      sums <- lapply(mcmcs,FUN=function(x){sum(x$X_sample[-(1:remove),t,j])})
      overall_mean[t,j] <- sum(as.numeric(as.character(sums)))/divisor
      acfs <- matrix(0,nrow=length(mcmcs),ncol=lag.max+1)
      for(i in 1:length(mcmcs)){
        setTxtProgressBar(pb, prog)
        prog <- prog + 1
        chain <- mcmcs[[i]]$X_sample[-(1:remove),t,j] - overall_mean[t,j]
        # print(length(chain))
        # print(length(acf(chain,demean=FALSE,lag.max=lag.max,type='covariance',plot=FALSE)$acf))
        acfs[i,] <- acf(chain,demean=FALSE,lag.max=lag.max,type='covariance',plot=FALSE)$acf
        
      }
      # average the autocovariance here
      average_acfs <- apply(acfs,MARGIN=2,FUN=mean)
      rho <- average_acfs/average_acfs[1]
      ac_time[t,j] <- (1+2*sum(rho[-1]))*tps
    }
  }
  close(pb)
  # return(ac_time)
  return(list(ac_time=ac_time,overall_mean=overall_mean))
}
