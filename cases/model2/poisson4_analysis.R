rm(list=ls())
load('./data/model2_poisson4.RData')
load('./data/model2_poisson4_env.RData')
load('./data/poisson4_lambda.RData')
load('./data/poisson4_pgmet.RData')
load('./data/poisson4_combined.RData')

start <- 200
N.plot <- 1000
N.lambda <- dim(lambda_out$X_sample)[1]
N.pgmet <- dim(pgmet_out$X_sample)[1]
N.combined <- dim(combined_out$X_sample)[1]
# index_lambda <- seq(start,N.lambda,round((N.lambda-start)/N.plot))
# index_pgmet <- seq(start,N.pgmet,round((N.pgmet-start)/N.plot))
# index_combined <- seq(start,N.combined,round((N.combined-start)/N.plot))
index_lambda <- sort(sample(N.lambda,N.plot))
index_pgmet <- sort(sample(N.pgmet,N.plot))
index_combined <- sort(sample(N.combined,N.plot))

cex.axis <- 2
cex.main <- 2
cex.lab <- 2
pch <- 16
cex <- 0.8

##### single state #####
t <- 250
j <- 4
pdf('./plots/poisson4_lambda_traceplot.pdf',width=11,height=8.5)
plot(lambda_out$X_sample[index_lambda,t,j],
     # main=expression('Trace plot of Lambda Sampler for X'['250,4']),
     ylab='',xlab='',ylim=c(-6,6),
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,cex=cex,pch=pch)
dev.off()
pdf('./plots/poisson4_pgmet_traceplot.pdf',width=11,height=8.5)
plot(pgmet_out$X_sample[index_pgmet,t,j],
     # main=expression('Trace plot of PGBS+Metropolis Sampler for X'['250,4']),
     ylab='',xlab='',ylim=c(-6,6),
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,cex=cex,pch=pch)
dev.off()
pdf('./plots/poisson4_combined_traceplot.pdf',width=11,height=8.5)
plot(combined_out$X_sample[index_combined,t,j],
     # main=expression('Trace plot of Combined Sampler for X'['250,4']),
     ylab='',xlab='',ylim=c(-6,6),
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,cex=cex,pch=pch)
dev.off()

##### state product ##### 
t1 <- 300
j1 <- 3
t2 <- 300
j2 <- 8

pdf('./plots/poisson4_lambda_prod_traceplot.pdf',width=11,height=8.5)
plot(lambda_out$X_sample[index_lambda,t1,j1]*lambda_out$X_sample[index_lambda,t2,j2],
     #main=expression('Trace plot of Lambda Sampler for X'['300,3']*'X'['300,8']),
     ylab='',xlab='',ylim=c(-15,20),
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,cex=cex,pch=pch)
dev.off()

pdf('./plots/poisson4_pgmet_prod_traceplot.pdf',width=11,height=8.5)
plot(pgmet_out$X_sample[index_pgmet,t1,j1]*pgmet_out$X_sample[index_pgmet,t2,j2],
     #main=expression('Trace plot of PGBS+Metropolis Sampler for X'['300,3']*'X'['300,8']),
     ylab='',xlab='',ylim=c(-15,20),
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,cex=cex,pch=pch)
dev.off()

pdf('./plots/poisson4_combined_prod_traceplot.pdf',width=11,height=8.5)
plot(combined_out$X_sample[index_combined,t1,j1]*combined_out$X_sample[index_combined,t2,j2],
     # main=expression('Trace plot of Combined Sampler for X'['300,3']*'X'['300,8']),
     ylab='',xlab='',ylim=c(-15,20),
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,cex=cex,pch=pch)
dev.off()


#### product state density comparison #### 
d1 <- density(lambda_out$X_sample[index_lambda,t1,j1]*lambda_out$X_sample[index_lambda,t2,j2])
d2 <- density(pgmet_out$X_sample[,t1,j1]*pgmet_out$X_sample[,t2,j2])
d3 <- density(combined_out$X_sample[index_combined,t1,j1]*combined_out$X_sample[index_combined,t2,j2])

# par(oma=c(2,10,0,4))
pdf('./plots/poisson4_prod_density.pdf',width=11,height=8.5)
par(mar=c(5.1,6.1,4.1,2.1))
plot(d1,ylim=c(0,0.11),main=expression('Density estimate for X'['300,3']*'X'['300,8']),
     cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,cex=cex,pch=pch)
lines(d2,col='red')
lines(d3,col='blue')
legend(9,0.10,legend=c('Lambda Sampler','PGMET Sampler','Combined Sampler'),
       col=c('black','red','blue'),cex=1.5,lty=1)
dev.off()

