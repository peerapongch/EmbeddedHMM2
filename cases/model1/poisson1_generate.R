# model 1 specification and simulation
seed <- 1
T <- 250; dim <- 3
mu_init <- rep(0,dim); rho <- 0.7; phi <- 0.9; phis <- rep(phi,dim); v <- 1/sqrt((1-phis^2)); sigma_init <- v %*% t(v)
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
delta <- rep(0.6,dim); c <- rep(-0.4,dim); F <- diag(phis); G <- t(chol(sigma)); Q <- diag(1,dim)
save.image(file='./data/model1_poisson1_env.RData')

source('../../models/model1_generate.R')
ssm_poisson <- generatePoissonGaussianSSM(T,dim,mu_init,sigma_init,F,G,Q,c,delta,seed)
save(ssm_poisson,file='./data/model1_poisson1.RData')