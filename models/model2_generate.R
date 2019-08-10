generate_gaussian_latent <- function(T,dim,mu_init,sigma_init,F,G,Q){
  require(MASS)
  # init
  sigma <- G %*% Q %*% t(G)
  # begin
  # x0 <- mvrnorm(1,mu_init,sigma_init)
  X <- matrix(0,ncol=dim,nrow=T)
  X[1,] <-  mvrnorm(1,mu_init,sigma_init)
  for(t in 2:T){
    X[t,] <- mvrnorm(1, F %*% X[t-1,], sigma) 
  }
  # return(list(x0=x0,X=X))
  return(X)
}

generate_model2_observation <- function(T,dim,X,delta){
  # as in model 1 in the paper 
  Y <- matrix(0,ncol=dim,nrow=T)
  abs_X <- abs(X)
  lambda <- t(apply(abs_X,MARGIN=1,FUN=function(x){x*delta}))
  for(t in 1:T){
    Y[t,] <- rpois(dim,lambda[t,])
  }
  return(Y)
}

generate_model2 <- function(T,dim,mu_init,sigma_init,F,G,Q,delta,seed){
  # call the above two functions 
  set.seed(seed)
  X <- generate_gaussian_latent(T,dim,mu_init,sigma_init,F,G,Q)
  Y <- generate_model2_observation(T,dim,X,delta)
  sigma <- G %*% Q %*% t(G)
  sigma_U <- chol(sigma)
  sigma_L <- t(sigma_U)
  sigma_inv <- chol2inv(sigma_U)
  sigma_init_U <- chol(sigma_init)
  sigma_init_L <- t(sigma_init_U)
  sigma_init_inv <- chol2inv(sigma_init_U)
  return(list(dim=dim,T=T,X=X,Y=Y,mu_init=mu_init,sigma_init=sigma_init,
              F=F,G=G,Q=Q,delta=delta,sigma=sigma,sigma_U=sigma_U,sigma_L=sigma_L,
              sigma_inv=sigma_inv,sigma_init_U=sigma_init_U,sigma_init_L=sigma_init_L,
              sigma_init_inv=sigma_init_inv, seed=seed))
}

