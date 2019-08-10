plot_mcmc_mu <- function(d,ssm,mcmc,mci_mu,plot.samples=FALSE,interval=10){
  # N <- mcmc$N
  N <- dim(mcmc$X_sample)[1]
  X_mcmc <- mcmc$X_sample
  plot(mci_mu[,d],col='black',type='l',ylim=c(-10,10),ylab=paste('X_',d,sep=''),xlab='t',
       main='Comparison of posterior smoothing means')
  if(plot.samples){
    for(i in seq(1,N,interval)){
      lines(X_mcmc[i,,d],col='grey')
    }
    lines(mci_mu[,d],col='black')
  }
  lines(ssm$X[,d],type='l',col='blue')
  legend(0,10,legend=c('Approximated smoothing mean','Latent state (truth)'),
         col=c('black','blue'), lty=c(1,1),cex=0.8)
}

plot_compare_mu <- function(mus,ssm,plot_index,ylim_upper=10,ylim_lower=-7,...){
  plot(ssm$X[,plot_index],type='l',ylim=c(ylim_lower,ylim_upper),...)
  lines(mus[[1]][1,,plot_index],col='red')
  lines(mus[[2]][1,,plot_index],col='blue')
  lines(mus[[3]][1,,plot_index],col='green')
  legend(0,ylim_upper,legend=c('True latent state','Metropolis','PGBS','PGBS+Metropolis'),
         col=c('black','red','blue','green'), lty=c(1,1),cex=0.8)
}
