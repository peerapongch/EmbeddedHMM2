setMKLthreads(2)
load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')

#### Random walk proposal ####
source('../../samplers/model1_lambda_param_RW.R')
N <- 9000
L <- 80
N.mcmc.param <- 20
seed <- 1 
rw_scale <- c(0.01,0.01,0.01,0.01)
checkpoint.name <- './data/checkpoint_poisson2_RW.RData'
filename <- './data/poisson2_lambda_param_RW.RData'
print(filename)
system.time(
  lambda_out <- lambdaModel1_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,
                                   seed=seed,rw_scale=rw_scale,checkpoint.name=checkpoint.name
  )
)

save(lambda_out,file=filename)

#### Random walk proposal with shift and scale update (RWSHSC1) #####
source('../../samplers/model1_lambda_param_RW_SHSC.R')
N <- 9000
L <- 80
N.mcmc.param <- 40
seed <- 1 
rw_scale <- c(0.01,0.01,0.01,0.01)
shsc_scale <- c(0.01,0.1)
filename <- './data/poisson2_lambda_param_RW_SHSC.RData'
print(filename)
checkpoint.name <- './data/checkpoint_RW_SHSC.RData'
system.time(
  lambda_out <- lambdaModel1_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,seed=seed,
                                   rw_scale=rw_scale,shsc_scale=shsc_scale,
                                   checkpoint.name=checkpoint.name
  )
)
save(lambda_out,file=filename)

#### Random walk proposal with shift and scale update without random walk update of sigma and c(RWSHSC2) #####
source('../../samplers/model1_lambda_param_RW_SHSC_2.R')
N <- 9000
L <- 80
N.mcmc.param <- 40
seed <- 1 
rw_scale <- c(0.01,0.01,0.01,0.01)
shsc_scale <- c(0.01,0.1)
filename <- './data/poisson2_lambda_param_RW_SHSC_2.RData'
print(filename)
checkpoint.name <- './data/checkpoint_RW_SHSC_2.RData'
system.time(
  lambda_out <- lambdaModel1_param(ssm_poisson,N,L,N.mcmc.param=N.mcmc.param,seed=seed,
                                   rw_scale=rw_scale,shsc_scale=shsc_scale,
                                   checkpoint.name=checkpoint.name
  )
)
save(lambda_out,file=filename)