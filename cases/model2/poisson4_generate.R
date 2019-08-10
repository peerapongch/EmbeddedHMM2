rm(list=ls())
# model 2 specification and simulation: poisson
seed <- 4
T <- 500; dim <- 10
mu_init <- rep(0,dim)
phi <- 0.9
rho <- 0.7
phis <- rep(phi,dim)
v <- 1/sqrt((1-phis^2))
sigma_init <- v %*% t(v)

for(i in 1:dim){
  for(j in 1:dim){
    if(i!=j){
      sigma_init[i,j] <- sigma_init[i,j] * rho
    }
  }
}
sigma <- matrix(rho,ncol=dim,nrow=dim)
for(i in 1:dim){
  sigma[i,i] <- 1
}

delta <- rep(0.8,dim)
F <- diag(phis)
G <- t(chol(sigma))
Q <- diag(1,dim)
save.image(file='./data/model2_poisson4_seed_env.RData')

source('../../models/model2_generate.R')
ssm_poisson <- generate_model2(T,dim,mu_init,sigma_init,F,G,Q,delta,seed)
# plot(ssm_poisson$Y[,1],type='l')
save(ssm_poisson,file='./data/model2_poisson4.RData')

# plot the observations
plot(ssm_poisson$Y[,1],type='l')
plot(ssm_poisson$Y[,2],type='l')
plot(ssm_poisson$Y[,3],type='l')
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
