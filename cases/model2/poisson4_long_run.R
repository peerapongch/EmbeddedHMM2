rm(list=ls())
setMKLthreads(2)
load('./data/model2_poisson4.RData')
load('./data/model2_poisson4_env.RData')

########## PGMET ##########
source('../../samplers/model2_pgmet.R')
N <- 550
es <- c(0.3,1)
L <- 5000
N.mcmc <- 50
seed <- 1
system.time(
  pgmet_out <- pgmetModel2(ssm_poisson,N,L,es,N.mcmc=N.mcmc,seed=seed)
)
save(pgmet_out,file='./data/poisson4_pgmet.RData')

########## Lambda ##########
source('../../samplers/model2_lambda.R')
N <- 36000
L <- 80
init <- matrix(1,nrow=T,ncol=dim)
system.time(
  lambda_out <- lambdaModel2(ssm_poisson,N,L,init=init)
)
save(lambda_out,file='./data/poisson4_lambda.RData')

########## Combined ##########
source('../../samplers/model2_combined.R')
N <- 500
L <- 80
L_particles <- 5000
es <- c(0.3,1)
init <- matrix(1,nrow=T,ncol=dim)
N.mcmc <- 50
system.time(
  combined_out <- combinedModel2(ssm_poisson,N,L,L_particles,es,N.mcmc=N.mcmc,init=init)
)
save(combined_out,file='./data/poisson4_combined.RData')