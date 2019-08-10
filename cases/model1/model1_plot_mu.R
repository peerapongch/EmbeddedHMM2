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
cex.main <- 1.2
cex.axis <- 1
cex.lab <- 1.2
colors <- c('black',brewer.pal(5,"Dark2"))
for(j in 1:dim){
  # jpeg(paste('poisson2_mu_x',j,'.jpg',sep=''),width = 750,height=550,quality=100)
  filename <- paste("./plots/poisson2_mu_dim",j,".pdf")
  pdf(file=filename,width=11,height=8.5)
  name <- paste('Smoothing mean of dimension ',j,sep='')
  plot(ssm_poisson$X[,j],type='l',main=name,ylim=c(-7,12),ylab='',xlab='Time',
       cex.main=1.5,cex.axis=1,cex.lab=1.2
       )
  lines(met_mu_avg[,j],col=colors[2])
  lines(pgbs_mu_avg[,j],col=colors[3])
  lines(pgmet_mu_avg[,j],col=colors[4])
  lines(lambda_mu_avg[,j],col=colors[5])
  lines(gamma_mu_avg[,j],col=colors[6])
  legend(0,12,legend = c("True latent state","Metropolis Sampler (200,000 samples)","PGBS Sampler (140,000 samples)",
                         "PGBS+Metropolis Sampler (280,000 samples)","Lambda Sampler (18,000 samples)",
                         "Gamma Sampler (40,000 samples)"),
         col=colors,cex=1.1,pt.cex=1,lty=1)
  dev.off()
}

