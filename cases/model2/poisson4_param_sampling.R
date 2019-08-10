# experiment with parameter sampling using different methods
rm(list=ls())
setMKLthreads(2)
load('./data/model2_poisson4.RData')
load('./data/model2_poisson4_env.RData')

#### Random walk metropolis ##### 
source('../../samplers/model2_lambda_param_RW.R')
N <- 9000
L <- 80
N.mcmc.param <- 20
seed <- 1 
rw_scale <- c(0.01,0.01,0.01)
checkpoint.name <- './data/checkpoint_poisson4_RW.RData'
filename <- './data/poisson4_lambda_param_RW.RData'
print(filename)
system.time(
  lambda_out <- lambdaModel2_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,
                                   seed=seed,rw_scale=rw_scale,checkpoint.name=checkpoint.name
  )
)

save(lambda_out,file=filename)

#### Random walk metropolis with scale update (RWSC1) #####
source('../../samplers/model2_lambda_param_RW_SC.R')
N <- 9000
L <- 80
N.mcmc.param <- 40
seed <- 1 
rw_scale <- c(0.01,0.01,0.01)
sc_scale <- 0.01
filename <- './data/poisson4_lambda_param_RW_SC.RData'
print(filename)
checkpoint.name <- './data/checkpoint_RW_SC.RData'
system.time(
  lambda_out <- lambdaModel2_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,seed=seed,
                                   rw_scale=rw_scale,sc_scale=sc_scale,
                                   checkpoint.name=checkpoint.name
  )
)

save(lambda_out,file=filename)

#### Random walk metropolis with scale update (RWSC2) with RW update of sigma #####
source('../../samplers/model2_lambda_param_RW_SC_2.R')
N <- 9000
L <- 80
N.mcmc.param <- 40
seed <- 1 
rw_scale <- c(0.01,0.01,0.01)
sc_scale <- 0.01
filename <- './data/poisson4_lambda_param_RW_SC_2.RData'
print(filename)
checkpoint.name <- './data/checkpoint_RW_SC_2.RData'
system.time(
  lambda_out <- lambdaModel2_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,seed=seed,
                                   rw_scale=rw_scale,sc_scale=sc_scale,
                                   checkpoint.name=checkpoint.name
  )
)

save(lambda_out,file=filename)