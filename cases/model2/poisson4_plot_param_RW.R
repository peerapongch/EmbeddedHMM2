library(coda)
load('./data/model2_poisson4.RData')
load('./data/model2_poisson4_env.RData')

load('./data/poisson4_lambda_param_RW.RData')

N <- dim(lambda_out$param_sample)[1]
start <- 10000
N.plot <- 1000
by <- round((N-start)/N.plot)
par(mar=c(5.1, 5.1, 4.1, 2.1))
cex.main <- 1.5
cex.axis <- 1.2
cex.lab <- 1.2
cex <- 0.7
pch <- 16
index <- seq(start,N,by)

pdf('./plots/poisson4_lambda_param_RW_phi.pdf',width=11,height=8.5)
plot(lambda_out$param_sample[index,1],xlab='Time',
     cex=cex,cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis,pch=pch,
     ylab=bquote(phi))
dev.off()

pdf('./plots/poisson4_lambda_param_RW_rho.pdf',width=11,height=8.5)
plot(lambda_out$param_sample[index,2],xlab='Time',
     cex=cex,cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis,pch=pch,
     ylab=bquote(rho))
dev.off()

pdf('./plots/poisson4_lambda_param_RW_sigma.pdf',width=11,height=8.5)
plot(lambda_out$param_sample[index,3],xlab='Time',
     cex=cex,cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis,pch=pch,
     ylab=bquote(sigma))
dev.off()