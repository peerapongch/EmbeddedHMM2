library(RColorBrewer)
load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
colors <- brewer.pal(dim,"Paired")

load('./data/poisson2_met_uactime.RData')
load('./data/poisson2_pgbs_uactime.RData')
load('./data/poisson2_pgmet_uactime.RData')
load('./data/poisson2_lambda_uactime.RData')
load('./data/poisson2_gamma_uactime.RData')
ylim_min <- 0
ylim_max <- 12.5

cex.main <- 1.5
cex.axis <- 2
cex.lab <- 2

# title <- 'Metropolis Sampler Adjusted Autocorrelation Time (0.05 s per sample)'
tps <- 0.05 
pdf(file="./plots/met_actime.pdf",width=11,height=8.5)
plot(met_uactime$ac_time[,1]*tps,type='l',col=colors[1],xlab='Time',cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='',ylim=c(ylim_min,ylim_max))
for(j in 2:dim){
  lines(met_uactime$ac_time[,j]*tps,col=colors[j])
}
abline(h=1,lty=2)
dev.off()

# title <- 'PGBS Sampler Adjusted Autocorrelation Time (0.19 s per sample)'
tps <- 0.19
pdf(file="./plots/pgbs_actime.pdf",width=11,height=8.5)
plot(pgbs_uactime$ac_time[,1]*tps,type='l',col=colors[1],xlab='Time',cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='',ylim=c(ylim_min,ylim_max))
for(j in 2:dim){
  lines(pgbs_uactime$ac_time[,j]*tps,col=colors[j])
}
abline(h=1,lty=2)
dev.off()


# title <- 'PGBS+Metropolis Sampler Adjusted Autocorrelation Time (0.12 s per sample)'
tps <- 0.12
pdf(file="./plots/pgmet_actime.pdf",width=11,height=8.5)
plot(pgmet_uactime$ac_time[,1]*tps,type='l',col=colors[1],xlab='Time',cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='',ylim=c(ylim_min,ylim_max))
for(j in 2:dim){
  lines(pgmet_uactime$ac_time[,j]*tps,col=colors[j])
}
abline(h=1,lty=2)
dev.off()

# title <- 'Lambda Sampler Adjusted Autocorrelation Time (0.54 s per sample)'
tps <- 0.54
pdf(file="./plots/lambda_actime.pdf",width=11,height=8.5)
plot(lambda_uactime$ac_time[,1]*tps,type='l',col=colors[1],xlab='Time',cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='',ylim=c(ylim_min,ylim_max))
for(j in 2:dim){
  lines(lambda_uactime$ac_time[,j]*tps,col=colors[j])
}
abline(h=1,lty=2)
dev.off()

# title <- 'Gamma Sampler Adjusted Autocorrelation Time (0.39 s per sample)'
ylim_min <- 0
ylim_max <- 36
tps <- 0.39
pdf(file="./plots/gamma_actime.pdf",width=11,height=8.5)
plot(gamma_uactime$ac_time[,1]*tps,type='l',col=colors[1],xlab='Time',cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='',ylim=c(ylim_min,ylim_max))
for(j in 2:dim){
  lines(gamma_uactime$ac_time[,j]*tps,col=colors[j])
}
abline(h=1,lty=2)
dev.off()

