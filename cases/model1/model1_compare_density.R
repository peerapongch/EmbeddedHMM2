load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')
library(coda)

source('../../samplers/model1_param_RW.R')
# N <- 50000
N <- 360000
rw_scale <- c(0.01,0.01,0.01,0.01)
seed <- 1
checkpoint.name <- './poisson1_param_RW.RData'
load('./data/poisson1_lambda_param_RW.RData')
X.init <- lambda_out$X_sample[lambda_out$N*2,,]
system.time(
  param_out <- paramModel1(ssm_poisson,N,X=X.init,rw_scale=rw_scale)
)

save(param_out,file='./data/poisson1_param_out_for_density_comparison.RData')

# if using presampled parameters
# load('./data/poisson1_lambda_param_RW.RData')
# load('./data/poisson1_param_out_for_density_comparison.RData')
cex.axis <- 2
cex.main <- 2
cex.lab <- 2
pch <- 16
cex <- 0.8

# compare density
pdf('./plots/poisson1_density_comparison_phi.pdf',width=11,height=8.5)
par(mar=c(5.1,5.1,4.1,2.1))
dense_cond <- density(param_out$param_sample[,1])
dense_marginal <- density(lambda_out$param_sample[,1])
plot(dense_marginal,ylim=c(0,max(dense_cond$y)),xlim=c(0.7,1),
     main='',
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab)
lines(dense_cond,col='red')
dev.off()

pdf('./plots/poisson1_density_comparison_rho.pdf',width=11,height=8.5)
par(mar=c(5.1,5.1,4.1,2.1))
dense_cond <- density(param_out$param_sample[,2])
dense_marginal <- density(lambda_out$param_sample[,2])
plot(dense_marginal,ylim=c(0,max(dense_cond$y)),xlim=c(0.2,0.9),
     main='',
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab)
lines(dense_cond,col='red')
dev.off()

pdf('./plots/poisson1_density_comparison_sigma.pdf',width=11,height=8.5)
par(mar=c(5.1,5.1,4.1,2.1))
dense_cond <- density(param_out$param_sample[,3])
dense_marginal <- density(lambda_out$param_sample[,3])
plot(dense_marginal,ylim=c(0,max(dense_cond$y)),xlim=c(0.3,0.9),
     main='',
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab)
lines(dense_cond,col='red')
dev.off()

pdf('./plots/poisson1_density_comparison_c.pdf',width=11,height=8.5)
par(mar=c(5.1,5.1,4.1,2.1))
dense_cond <- density(param_out$param_sample[,4])
dense_marginal <- density(lambda_out$param_sample[,4])
plot(dense_marginal,ylim=c(0,max(dense_cond$y)),xlim=c(-1.1,1),
     main='',
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab)
lines(dense_cond,col='red')
dev.off()
