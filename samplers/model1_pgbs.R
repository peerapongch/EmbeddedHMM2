forward_step <- function(X_current,Y,T,L,dim,F,sigma_U,mu_init,sigma_init,c,delta){
  X_pool <- array(0,dim=c(T,L,dim))
  W <- matrix(0,nrow=T,ncol=L)
  lW <- matrix(0,nrow=T,ncol=L)
  for(t in 1:T){
    X_pool[t,1,] <- X_current[t,]
    if(t==1){
      X_pool[t,2:L,] <- mvrnorm(L-1,mu_init,sigma_init)
    } else {
      anc <- sample(1:L,L-1,replace=TRUE,prob=W[t-1,])
      Z <- matrix(rnorm(dim*(L-1)),ncol=dim,nrow=L-1)
      X_pool[t,2:L,] <- X_pool[t-1,anc,]%*%t(F) + Z %*% sigma_U
    }
    lambdas <- t(X_pool[t,,])*delta+c
    lW[t,] <- colSums(dpois(Y[t,],exp(lambdas),log=TRUE))
    num <- exp(lW[t,]-max(lW[t,]))
    W[t,] <- num/sum(num)
  }
  return(list(lW=lW,W=W,X_pool=X_pool))
}

backward_step <- function(forward_results,T,L,F,sigma_inv){
  W <- forward_results$W 
  lW <- forward_results$lW 
  X_pool <- forward_results$X_pool
  X_new <- matrix(0,ncol=dim,nrow=T)
  l_T <- sample(1:L,1,replace=TRUE,prob=W[T,])
  X_new[T,] <- X_pool[T,l_T,]
  for(t in (T-1):1){
    diff <- X_new[t+1,] - t(X_pool[t,,] %*% F) # also only for symmetric F 
    lprob <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    lprob <- lprob+lW[t,]
    num <- exp(lprob-max(lprob))
    prob <- num/sum(num)
    # setting new x
    X_new[t,] <- X_pool[t,sample(1:L,1,prob=prob),]
  } 
  return(X_new)
}

pgbsModel1 <- function(ssm,N,L,init=NULL,seed=NULL,return.weight=FALSE){
  require(MASS)
  #poisson observation and gaussian latent process 
  mu_init <- ssm$mu_init
  sigma_init <- ssm$sigma_init
  F <- ssm$F
  c <- ssm$c
  delta <- ssm$delta
  T <- ssm$T
  sigma_U <- ssm$sigma_U
  sigma_L <- ssm$sigma_L
  Y <- ssm$Y
  dim <- ssm$dim
  
  if(is.null(ssm$sigma_inv)){
    print('old ssm object, computing inverse for sigma')
    sigma_inv <- chol2inv(sigma_U)
  } else {
    sigma_inv <- ssm$sigma_inv
  }

  X_sample <- array(0,dim=c(2*N+1,T,dim))
  if(return.weight){
    W <- array(logical(0),dim=c(2*N,T,L))
    lW <- array(logical(0),dim=c(2*N,T,L))
  }

  if(is.null(init)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    X_sample[1,,] <- mvrnorm(T,mu_init,sigma_init)
  } else {
    X_sample[1,,] <- init 
  }
  X_current <- X_sample[1,,]
  
  pb <- txtProgressBar(min=0,max=N*2,title="pgbs",style=3)
  acceptance_rate <- 0
  for(n in seq(2,2*N,2)){
    ### forward sequence ###
    forward_results <- forward_step(X_current,Y,T,L,dim,F,sigma_U,mu_init,sigma_init,c,delta)
    if(return.weight){
      W[n-1,,] <- forward_results$W
      lW[n-1,,] <- forward_results$lW
    }
    backward_out <- backward_step(forward_results,T,L,F,sigma_inv)
    if(any(backward_out != X_current)){
      X_current <- backward_out
      acceptance_rate <- acceptance_rate + 1
    }
    # save 
    X_sample[n,,] <- X_current
    
    ### reversed sequence ###
    forward_results <- forward_step(X_current[seq(T,1,-1),],Y[seq(T,1,-1),],T,L,dim,F,sigma_U,mu_init,sigma_init,c,delta)
    if(return.weight){
      W[n,seq(T,1,-1),] <- forward_results$W
      lW[n,seq(T,1,-1),] <- forward_results$lW
    }
    backward_out <- backward_step(forward_results,T,L,F,sigma_inv)
    if(any(backward_out[seq(T,1,-1),] != X_current)){
      X_current <- backward_out[seq(T,1,-1),]
      acceptance_rate <- acceptance_rate + 1
    }
    # save
    X_sample[n+1,,] <- X_current
    
    setTxtProgressBar(pb, n)
  }
  close(pb)
  acceptance_rate <- acceptance_rate/(2*N)
  if(return.weight){
    return(list(X_sample = X_sample[-1,,],N=N,init=init,W=W,lW=lW,seed=seed,acceptance_rate=acceptance_rate))
  }
  return(list(X_sample = X_sample[-1,,],N=N,init=init,seed=seed,acceptance_rate=acceptance_rate))
}