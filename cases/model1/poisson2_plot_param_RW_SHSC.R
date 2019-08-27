rm(list=ls())
library(coda)
load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')

load('./data/poisson2_lambda_param_RW_SHSC.RData')

N <- dim(lambda_out$param_sample)[1]
start <- 10000
N.plot <- 1000
by <- round((N-start)/N.plot)
par(mar=c(5.1, 6.1, 4.1, 2.1))
cex.axis <- 2
cex.main <- 2
cex.lab <- 2
pch <- 16
cex <- 0.8
index <- seq(start,N,by)

pdf('./plots/poisson2_lambda_param_RW_SHSC_phi.pdf',width=11,height=8.5)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(lambda_out$param_sample[index,1],xlab='',
     cex=cex,cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis,pch=pch,
     ylab=bquote(phi),ylim=c(0.82,0.94))
dev.off()

pdf('./plots/poisson2_lambda_param_RW_SHSC_rho.pdf',width=11,height=8.5)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(lambda_out$param_sample[index,2],xlab='',
     cex=cex,cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis,pch=pch,
     ylab=bquote(rho),ylim=c(0.45,0.8))
dev.off()

pdf('./plots/poisson2_lambda_param_RW_SHSC_sigma.pdf',width=11,height=8.5)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(lambda_out$param_sample[index,3],xlab='',
     cex=cex,cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis,pch=pch,
     ylab=bquote(sigma),ylim=c(0.52,0.78))
dev.off()

pdf('./plots/poisson2_lambda_param_RW_SHSC_c.pdf',width=11,height=8.5)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(lambda_out$param_sample[index,4],xlab='',
     cex=cex,cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis,pch=pch,
     ylab=bquote(c),ylim=c(-2,0.5))
dev.off()
