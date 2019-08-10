shift_update <- function(this_X_current,next_X_pool,Y,t,l,l_current,rev_F,c,delta,Us,l_new){
  # shift update: update l then x
  # update of l UAR(1:L) done above since independent of l
  next_X_new <- next_X_pool[l_new[l],]
  next_X_current <- next_X_pool[l_current,]
  X_new <- this_X_current + rev_F %*% (next_X_new-next_X_current)
  lnum <- sum(dpois(Y[t+1,],exp(c+delta*next_X_new),log=TRUE))
  ldenom <- sum(dpois(Y[t+1,],exp(c+delta*next_X_current),log=TRUE))
  hastings <- exp(lnum - ldenom)
  alpha <- min(1,hastings)
  if(alpha>Us[2,l]){
    return(list(X_new=X_new,l_new=l_new[l]))
  } else {
    return(list(X_new=this_X_current,l_new=l_current))
  }
}

backward_pool <- function(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,c,delta,rev_F,rev_sigma_U,elim){
  # metropolis
  # sequential update
  X_pool <- array(0,dim=c(T,L,dim))
  acceptance_rate <- matrix(0,nrow=T,ncol=2)
  # start
  es <- runif(L,elim[1],elim[2])
  es_2 <- es^2
  zs <- matrix(rnorm(L*dim),ncol=dim,nrow=L)
  Us <- runif(L)
  
  # sample the index for the current states 
  k <- sample(1:L,T,replace=TRUE)
  
  # place the first current state 
  X_pool[T,k[T],] <- X_current[T,]
  
  # do for terminal T 
  # transition down from k[T] to 1
  if(k[T]>1){
    this_X_current <- X_current[T,]
    for(l in (k[T]-1):1){
      # update latent states
      # autoregressive update 
      this_X_current <- mu_init + sqrt(1-es_2[l])*(this_X_current-mu_init) + es[l]*sigma_init_L%*%zs[l,]
      X_pool[T,l,] <- this_X_current
    } 
  }
  
  if(k[T]<L){
    this_X_current <- X_current[T,]
    for(l in (k[T]+1):L){
      # autoregressive update latent states
      this_X_current <- mu_init + sqrt(1-es_2[l])*(this_X_current-mu_init) + es[l]*sigma_init_L%*%zs[l,]
      X_pool[T,l,] <- this_X_current
    }
  }
  acceptance_rate[T,] <- c(1,0)

  # then for t>1
  for(t in (T-1):1){
    this_X_current <- X_current[t,]
    #### stochastic initialisation of l 
    diff <- t(X_pool[t+1,,]) - as.vector(F %*% this_X_current)
    lprob_trans <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    lprob_obs <- colSums(dpois(Y[t+1,],exp(c+delta*t(X_pool[t+1,,])),log=TRUE))
    lprob <- lprob_obs + lprob_trans
    prob <- exp(lprob-max(lprob)); prob <- prob/sum(prob)
    l_original <- sample(1:L,1,prob=prob)
    ####

    this_X_pool <- matrix(logical(0),nrow=L,ncol=dim) # for this timestep only 
    this_X_pool[k[t],] <- this_X_current
    Us <- matrix(runif(L*2),nrow=2,ncol=L)
    l_new <- sample(1:L,L,replace=TRUE)

    # setup for autoregressive update
    mu <- X_pool[t+1,,] %*% rev_F # symmetric F  
    es <- runif(L,elim[1],elim[2]) # hardcoded
    sq_term <- sqrt(1-es^2)
    zs <- matrix(rnorm(L*dim),nrow=L,ncol=dim) 
    last <- es*zs%*%rev_sigma_U
    
    x_accept <- 0 # to measure acceptance rates 
    l_accept <- 0

    ### begin
    next_X_pool <- X_pool[t+1,,]
    # reversed transition
    l_current <- l_original
    if(k[t]>1){
      this_X_current <- X_current[t,]
      for(l in (k[t]-1):1){
        # Shift
        shift_out <- shift_update(this_X_current,next_X_pool,Y,t,l,l_current,rev_F,c,delta,Us,l_new)
        # count and then update
        if(all(this_X_current!=shift_out$X_new)){
          l_accept <- l_accept + 1
        }
        this_X_current <- shift_out$X_new
        l_current <- shift_out$l_new

        # AR
        X_new <- mu[l_current,] + sq_term[l]*(this_X_current-mu[l_current,]) + last[l,]
        if(all(this_X_current!=X_new)){
          x_accept <- x_accept + 1  
        }      
        this_X_current <- X_new
       
        # set pool
        this_X_pool[l,] <- this_X_current
      } 
    }

    # forward transition
    l_current <- l_original
    if(k[t]<L){
      this_X_current <- X_current[t,]
      for(l in (k[t]+1):L){
        # AR 
        X_new <- mu[l_current,] + sq_term[l]*(this_X_current-mu[l_current,]) + last[l,]

        if(all(this_X_current!=X_new)){
          x_accept <- x_accept + 1 
        }      
        this_X_current <- X_new
        
        # Shift
        shift_out <- shift_update(this_X_current,next_X_pool,Y,t,l,l_current,rev_F,c,delta,Us,l_new)
        # count and then update
        if(all(this_X_current!=shift_out$X_new)){
          l_accept <- l_accept + 1
        }
        this_X_current <- shift_out$X_new
        l_current <- shift_out$l_new

        # set pool
        this_X_pool[l,] <- this_X_current
      } 
    }
    
    acceptance_rate[t,] <- c(x_accept,l_accept)/(L-1)
    X_pool[t,,] <- this_X_pool
    X_current[t,] <- this_X_current

  } 
  return(list(X_pool=X_pool,acceptance_rate=acceptance_rate))
}

forward_sampling <- function(X_pool,Y,L,T,dim,F,mu_init,sigma_init_inv,sigma_inv){
  # backward sampling
  X_new <- matrix(0,nrow=T,ncol=dim)
  # 2 parts: prior and observation model 
  # prior 
  diff_1 <- t(X_pool[1,,]) - mu_init
  lprob_prior <- diag(-1/2*t(diff_1)%*%sigma_init_inv%*%diff_1) # expect length L
  # obs
  lprob_obs <- colSums(dpois(Y[1,],exp(c+delta*t(X_pool[1,,])),log=TRUE))
  # total
  lprob <- lprob_obs + lprob_prior
  prob <- exp(lprob-max(lprob)); prob <- prob/sum(prob)
  X_new[1,] <- X_pool[1,sample(1:L,1,prob=prob),]
  for(t in 2:T){
    diff <- t(X_pool[t,,]) - as.vector(F %*% X_new[t-1,] )
    lprob_trans <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    lprob_obs <- colSums(dpois(Y[t,],exp(c+delta*t(X_pool[t,,])),log=TRUE))
    lprob <- lprob_obs + lprob_trans
    prob <- exp(lprob-max(lprob)); prob <- prob/sum(prob)
    X_new[t,] <- X_pool[t,sample(1:L,1,prob=prob),]
  }
  return(X_new)
}

gammaModel1 <- function(ssm,N,L,init=NULL,seed=NULL,elim=c(0,0.4)){
  # sampling epsilon values instead
  # setup
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
  # for backward compatibility 
  if((is.null(ssm$sigma_init_L)) | (is.null(ssm$sigma_init_U)) | (is.null(ssm$sigma_init_inv))){
    print('old ssm object, computing Choleskey and inverse for sigma_init')
    sigma_init_U <- chol(sigma_init)
    sigma_init_L <- t(sigma_init_U)
    sigma_init_inv <- chol2inv(sigma_init_U)
  } else {
    sigma_init_L <- ssm$sigma_init_L
    sigma_init_U <- ssm$sigma_init_U
    sigma_init_inv <- ssm$sigma_init_inv
  }
  if(is.null(ssm$sigma_inv)){
    print('old ssm object, computing inverse for sigma')
    sigma_inv <- chol2inv(sigma_U)
  } else {
    sigma_inv <- ssm$sigma_inv
  }
  rev_F <- solve(F)
  rev_sigma <- rev_F %*% sigma %*% t(rev_F)
  rev_sigma_U <- chol(rev_sigma)

  ### begin
  X_sample <- array(0,dim=c(2*N+1,T,dim))
  if(is.null(init)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    X_sample[1,,] <- mvrnorm(T,mu_init,sigma_init)
  } else { 
    X_sample[1,,] <- init 
  }
  X_current <- X_sample[1,,]

  acceptance_rate <- array(0,dim=c(2*N,T,2))
  pb <- txtProgressBar(min=0,max=2*N,title="ehmm",style=3)

  for(i in seq(2,2*N,2)){
    # forward sequence
    pool_out <- backward_pool(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,c,delta,rev_F,rev_sigma_U,elim)
    X_pool <- pool_out$X_pool
    acceptance_rate[i-1,,] <- pool_out$acceptance_rate
    
    X_current <- forward_sampling(X_pool,Y,L,T,dim,F,mu_init,sigma_init_inv,sigma_inv)
    X_sample[i,,] <- X_current
    
    # reversed sequence
    pool_out <- backward_pool(X_current[seq(T,1,-1),],Y[seq(T,1,-1),],T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,c,delta,rev_F,rev_sigma_U,elim)
    X_pool <- pool_out$X_pool
    acceptance_rate[i,seq(T,1,-1),] <- pool_out$acceptance_rate

    X_current[seq(T,1,-1),] <- forward_sampling(X_pool,Y[seq(T,1,-1),],L,T,dim,F,mu_init,sigma_init_inv,sigma_inv)
    X_sample[i+1,,] <- X_current
    
    setTxtProgressBar(pb, i)
  }
  return(list(X_sample=X_sample[-1,,],N=N,L=L,init=init,seed=seed,X_pool=X_pool,
              acceptance_rate=acceptance_rate))
}