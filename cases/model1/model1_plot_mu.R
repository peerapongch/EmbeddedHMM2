load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
library(RColorBrewer)

load('./data/poisson2_met_mu.RData')
met_mu_avg <- matrix(logical(0),nrow=T,ncol=dim)
for(t in 1:T){
  met_mu_avg[t,] <- colSums(met_mu[,t,])/dim(met_mu[,t,])[1]
}

load('./data/poisson2_pgbs_mu.RData')
pgbs_mu_avg <- matrix(logical(0),nrow=T,ncol=dim)
for(t in 1:T){
  pgbs_mu_avg[t,] <- colSums(pgbs_mu[,t,])/dim(pgbs_mu[,t,])[1]
}

load('./data/poisson2_pgmet_mu.RData')
pgmet_mu_avg <- matrix(logical(0),nrow=T,ncol=dim)
for(t in 1:T){
  pgmet_mu_avg[t,] <- colSums(pgmet_mu[,t,])/dim(pgmet_mu[,t,])[1]
}

load('./data/poisson2_lambda_mu.RData')
lambda_mu_avg <- matrix(logical(0),nrow=T,ncol=dim)
for(t in 1:T){
  lambda_mu_avg[t,] <- colSums(lambda_mu[,t,])/dim(lambda_mu[,t,])[1]
}

load('./data/poisson2_gamma_mu.RData')
gamma_mu_avg <- matrix(logical(0),nrow=T,ncol=dim)
for(t in 1:T){
  gamma_mu_avg[t,] <- colSums(gamma_mu[,t,])/dim(gamma_mu[,t,])[1]
}
par(mar=c(5.1, 6.1, 4.1, 2.1))
cex.axis <- 2
cex.main <- 2
cex.lab <- 2
pch <- 16
cex <- 0.8
colors <- c('black',brewer.pal(5,"Dark2"))
for(j in 1:dim){
  filename <- paste("./plots/poisson2_mu_dim",j,".pdf",sep='')
  pdf(file=filename,width=11,height=8.5)
  # name <- paste('Smoothing mean of dimension ',j,sep='')
  plot(ssm_poisson$X[,j],type='l', #main=name,
       ylim=c(-7,14),ylab='',xlab='Time',
       cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab
       )
  lines(met_mu_avg[,j],col=colors[2])
  lines(pgbs_mu_avg[,j],col=colors[3])
  lines(pgmet_mu_avg[,j],col=colors[4])
  lines(lambda_mu_avg[,j],col=colors[5])
  lines(gamma_mu_avg[,j],col=colors[6])
  legend(0,14,legend = c("True latent state","Metropolis Sampler","PGBS Sampler",
                         "PGMET Sampler","Lambda Sampler",
                         "Gamma Sampler"),
         col=colors,cex=1.6,lty=1)
  dev.off()
}

