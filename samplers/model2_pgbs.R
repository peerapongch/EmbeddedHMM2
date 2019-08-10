forward_step <- function(X_current,Y,T,L,dim,F,sigma_U,mu_init,sigma_init,delta){
  X_pool <- array(0,dim=c(T,L,dim))
  W <- matrix(0,nrow=T,ncol=L)
  lW <- matrix(0,nrow=T,ncol=L)
  this_X_pool <- matrix(logical(0),nrow=L,ncol=dim)
  for(t in 1:T){
    # X_pool[t,1,] <- X_sample[n-1,t,]
    this_X_pool[1,] <- X_current[t,]
    # X_pool[t,1,] <- X_current[t,]
    if(t==1){
      this_X_pool[2:L,] <- mvrnorm(L-1,mu_init,sigma_init)
      # X_pool[t,2:L,] <- mvrnorm(L-1,mu_init,sigma_init)
    } else {
      anc <- sample(1:L,L-1,replace=TRUE,prob=W[t-1,])
      Z <- matrix(rnorm(dim*(L-1)),ncol=dim,nrow=L-1)
      this_X_pool[2:L,] <- X_pool[t-1,anc,]%*%t(F) + Z %*% sigma_U 
      # X_pool[t,2:L,] <- X_pool[t-1,anc,]%*%t(F) + Z %*% sigma_U
    }

    lambdas <- delta*abs(t(this_X_pool))
    lW[t,] <- colSums(dpois(Y[t,],lambdas,log=TRUE))
    num <- exp(lW[t,]-max(lW[t,]))
    W[t,] <- num/sum(num)
    
    # set pool
    X_pool[t,,] <- this_X_pool
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
    # diff <- X_new[t+1,] - X_pool[t,,] %*% F
    diff <- X_new[t+1,] - t(X_pool[t,,] %*% F) # also only for symmetric F 
    lprob <- diag(-1/2*t(diff)%*%sigma_inv%*%diff) # consider inverse once, consider debug (refer ehmm line 76): done
    # lprob <- apply(X_pool[t,,],MARGIN=1,FUN=ld_model1_transition,x_to=X_new[t+1,],F=F,sigma_U=sigma_U)
    lprob <- lprob+lW[t,]
    num <- exp(lprob-max(lprob))
    prob <- num/sum(num)
    # setting new x
    X_new[t,] <- X_pool[t,sample(1:L,1,prob=prob),]
  } 
  return(X_new)
}

pgbsModel2 <- function(ssm,N,L,init=NULL,seed=NULL){
  require(MASS)
  #poisson observation and gaussian latent process 
  mu_init <- ssm$mu_init
  sigma_init <- ssm$sigma_init
  F <- ssm$F
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
  if(is.null(init)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    X_sample[1,,] <- mvrnorm(T,mu_init,sigma_init)
  } else {
    if(!is.null(seed)){
      set.seed(seed)
    }
    X_sample[1,,] <- init 
  }
  X_current <- X_sample[1,,]
  
  pb <- txtProgressBar(min=0,max=N*2,title="pgbs",style=3)
  for(n in seq(2,2*N,2)){
    
    ### forward sequence ###
    forward_results <- forward_step(X_current,Y,T,L,dim,F,sigma_U,mu_init,sigma_init,delta)
    X_new <- backward_step(forward_results,T,L,F,sigma_inv)
    # save 
    X_sample[n,,] <- X_new
    
    ### reversed sequence ###
    X_current <- X_new[seq(T,1,-1),] # form reversed sequence
    forward_results <- forward_step(X_current,Y[seq(T,1,-1),],T,L,dim,F,sigma_U,mu_init,sigma_init,delta)
    X_new <- backward_step(forward_results,T,L,F,sigma_inv)
    # save
    X_sample[n+1,,] <- X_new[seq(T,1,-1),] # reverse the reversed
    
    ### update current and terminate ###
    X_current <- X_sample[n+1,,] # forward sequence
    setTxtProgressBar(pb, n)
  }
  close(pb)
  return(list(X_sample = X_sample[-1,,],N=N,init=init,
              seed=seed))
}