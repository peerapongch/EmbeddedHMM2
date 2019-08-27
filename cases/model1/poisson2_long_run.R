setMKLthreads(2)
load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')

# seeds for sampling
seeds <- 1:5

######### Metropolis #########
source('../../samplers/model1_met.R')
N <- 3000000
es <- c(0.2,0.8)
for(s in seeds){
  seed <- s
  filename <- paste('./data/poisson2_met_seed',seed,'.RData',sep='')
  print(filename)
  system.time(
    met_out <- metModel1(ssm_poisson,N,es,seed=seed) 
  )
  save(met_out,file=filename)
}

######### PGBS #########
source('../../samplers/model1_pgbs.R')
N <- 40000
L <- 250
for(s in seeds){
  seed <- s
  filename <- paste('./data/poisson2_pgbs_seed',seed,'.RData',sep='')
  print(filename)
  system.time(
    pgbs_out <- pgbsModel1(ssm_poisson,N,L,seed=seed) 
  )
  save(pgbs_out,file=filename)
}

######### PGBS+Metropolis #########
source('../../samplers/model1_pgmet.R')
N <- 30000
L <- 250
es <- c(0.2,0.8)
for(s in seeds){
  seed <- s
  filename <- paste('./data/poisson2_pgmet_seed',seed,'.RData',sep='')
  print(filename)
  system.time(
    pgmet_out <- pgmetModel1(ssm_poisson,N,L,es,seed=seed) 
  )
  save(pgmet_out,file=filename)
}

######### Lambda #########
source('../../samplers/model1_lambda.R')
N <- 14000
L <- 50
for(s in seeds){
  seed <- s
  filename <- paste('./data/poisson2_lambda_seed',seed,'.RData',sep='')
  print(filename)
  system.time(
    lambda_out <- lambdaModel1(ssm_poisson,N,L,seed=seed) 
  )
  save(lambda_out,file=filename)
}

######### Gamma #########
source('../../samplers/model1_gamma.R')
N <- 20000
L <- 50
elim <- c(0,0.2)
for(s in seeds){
  seed <- s
  filename <- paste('./data/poisson2_gamma_seed',seed,'.RData',sep='')
  print(filename)
  system.time(
    gamma_out <- gammaModel1(ssm_poisson,N,L,seed=seed,elim=elim) 
  )
  save(gamma_out,file=filename)
}

