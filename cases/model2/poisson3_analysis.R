rm(list=ls())
load('./data/model2_poisson3.RData')
load('./data/model2_poisson3_env.RData')
load('./data/poisson3_lambda.RData')
load('./data/poisson3_pgmet.RData')
load('./data/poisson3_combined.RData')

start <- 1000
N.plot <- 2000
N.lambda <- dim(lambda_out$X_sample)[1]
N.pgmet <- dim(pgmet_out$X_sample)[1]
N.combined <- dim(combined_out$X_sample)[1]
index_lambda <- seq(start,N.lambda,round((N.lambda-start)/N.plot))
index_pgmet <- seq(start,N.pgmet,round((N.pgmet-start)/N.plot))
index_combined <- seq(start,N.combined,round((N.combined-start)/N.plot))

### analyse trace plot, single
t <- 100
j <- 3
plot(lambda_out$X_sample[index_lambda,t,j],main=bquote('Trace plot of Lambda Sampler for X'['100,3']),ylab='',xlab='Time')
plot(pgmet_out$X_sample[index_pgmet,t,j],main=bquote('Trace plot of PGBS+Metropolis Sampler for X'['100,3']),ylab='',xlab='Time')
plot(combined_out$X_sample[index_combined,t,j],main=bquote('Trace plot of Combined Sampler for X'['100,3']),ylab='',xlab='Time')

### analyse trace plot, product
t1 <- 100
j1 <- 1
t2 <- 100
j2 <- 3
plot(lambda_out$X_sample[index_lambda,t1,j1],main='Trace plot of Lambda Sampler for X_100,1',ylab='')
plot(lambda_out$X_sample[index_lambda,t2,j2],main='Trace plot of Lambda Sampler for X_100,3',ylab='')
plot(lambda_out$X_sample[index_lambda,t1,j1]*lambda_out$X_sample[index_lambda,t2,j2],main='Trace plot of product of states, Lambda Sampler',ylab='')

plot(pgmet_out$X_sample[index_pgmet,t1,j1],main='Trace plot of PGBS+Metropolis Sampler for X_100,1',ylab='')
plot(pgmet_out$X_sample[index_pgmet,t2,j2],main='Trace plot of PGBS+Metropolis Sampler for X_100,3',ylab='')
plot(pgmet_out$X_sample[index_pgmet,t1,j1]*pgmet_out$X_sample[index_pgmet,t2,j2],main='Trace plot of product of states, PGBS+Metropolis',ylab='')

plot(combined_out$X_sample[index_combined,t1,j1],main='Trace plot of Combined Sampler for X_100,1',ylab='')
plot(combined_out$X_sample[index_combined,t2,j2],main='Trace plot of Combined Sampler for X_100,3',ylab='')
plot(combined_out$X_sample[index_combined,t1,j1]*combined_out$X_sample[index_combined,t2,j2],main='Trace plot of product of states, Combined Sampler',ylab='')

### positivity
t <- 100
plot(pgmet_out$X_sample[,t,1],pgmet_out$X_sample[,t,2],xlab='X_100,1',ylab='x_100,2')
plot(pgmet_out$X_sample[,t,1],pgmet_out$X_sample[,t,3],xlab='X_100,1',ylab='x_100,3')
plot(pgmet_out$X_sample[,t,2],pgmet_out$X_sample[,t,3],xlab='X_100,2',ylab='x_100,3')
pairs(pgmet_out$X_sample[,100,])

### analyse smoothing mean
source('../../evaluation/function_mcmc_mu.R')
lambda_mu <- mci_mu(lambda_out$X_sample,T,dim)
pgmet_mu <- mci_mu(pgmet_out$X_sample,T,dim)
combined_mu <- mci_mu(combined_out$X_sample,T,dim)

library(RColorBrewer)
colors <- c('black',brewer.pal(3,"Dark2"))
for(j in 1:dim){
  filename <- paste("./plots/poisson3_mu_dim",j,'.pdf',sep='')
  pdf(file=filename,width=11,height=8.5)
  title <- paste('Smoothing mean of dimension ',j,sep='')
  plot(ssm_poisson$X[,j],type='l',col=colors[1],ylim=c(-8,12),main=title,ylab='',xlab='Time',
       cex.main=1.5,cex.axis=1,cex.lab=1.2
       )
  lines(lambda_mu[,j],col=colors[2])
  lines(pgmet_mu[,j],col=colors[3])
  lines(combined_mu[,j],col=colors[4])
  legend(300,12,legend=c("True value","Lambda Sampler","PGBS+Metropolis Sampler","Combined Sampler"),
         cex=1,lty=1,col=colors)
  dev.off()
}

### more traceplots 
library(coda)
load('./data/poisson3_pgmet_v4.RData')
x1 <- as.mcmc(pgmet_out$X_sample[,1,])
load('./data/poisson3_pgmet_v2.RData')
x2 <- as.mcmc(pgmet_out$X_sample[1:500,100,])
plot(x1)
plot(x2)
load('./data/poisson3_lambda.RData')
x3 <- as.mcmc(lambda_out$X_sample[,1,])
plot(x3)
# x <- as.mcmc(pgmet_out$X_sample[,100,])
# plot(x)
# x <- as.mcmc(lambda_out$X_sample[1:5000,100,])
# plot(x)
