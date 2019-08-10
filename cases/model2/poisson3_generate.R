rm(list=ls())
# model 2 specification and simulation
seed <- 4
T <- 500; dim <- 3
mu_init <- rep(0,dim)
phi <- 0.9
rho <- 0.7
phis <- rep(phi,dim)
v <- 1/sqrt((1-phis^2))
sigma_init <- rho * v %*% t(v)
diag(sigma_init) <- diag(sigma_init)/rho

sigma <- matrix(rho,ncol=dim,nrow=dim)
for(i in 1:dim){
  sigma[i,i] <- 1
}

delta <- rep(0.8,dim)
F <- diag(phis)
G <- t(chol(sigma))
Q <- diag(1,dim)
save.image(file='./data/model2_poisson3_env.RData')

source('../../models/model2_generate.R')
ssm_poisson <- generate_model2(T,dim,mu_init,sigma_init,F,G,Q,delta,seed)
save(ssm_poisson,file='./data/model2_poisson3.RData')

# plot the observations
plot(ssm_poisson$Y[,1],type='l')
plot(ssm_poisson$Y[,2],type='l')
plot(ssm_poisson$Y[,3],type='l')
apply(ssm_poisson$Y,MARGIN=2,FUN=mean)

load('./data/model2_poisson3_seed3.RData')
plot(ssm_poisson$Y[,1],type='l')
plot(ssm_poisson$Y[,2],type='l')
plot(ssm_poisson$Y[,3],type='l')
apply(ssm_poisson$Y,MARGIN=2,FUN=mean)
# plot(ssm_poisson$Y[,4],type='l')
# plot(ssm_poisson$Y[,5],type='l')

plot(ssm_poisson$X[,1],type='l')
lines(ssm_poisson$X[,2])
lines(ssm_poisson$X[,3],col='red')
# lines(ssm_poisson$X[,4])
# lines(ssm_poisson$X[,5])

# plot the lambdas upto delta
plot(abs(ssm_poisson$X[,1]),type='l')
lines(abs(ssm_poisson$X[,2]))
lines(abs(ssm_poisson$X[,3]),col='red')
# lines(abs(ssm_poisson$X[,4]))
# lines(abs(ssm_poisson$X[,5]))

plot(ssm_poisson$Y[,1],type='l')
lines(abs(ssm_poisson$X[,1]),col='red')
