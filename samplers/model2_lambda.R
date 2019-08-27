autoregressive_update <- function(this_X_current,mu,Y,t,l,l_current,F,delta,sq_term,Us,last){
  mu_current <- mu[l_current,]
  X_new <- mu_current + sq_term[l]*(this_X_current-mu_current)+last[l,]
  # hastings ratio
  lnum <- sum(dpois(Y[t,],delta*abs(X_new),log=TRUE))
  ldenom <- sum(dpois(Y[t,],delta*abs(this_X_current),log=TRUE))
  lhastings <- lnum - ldenom
  hastings <- exp(lhastings)
  alpha <- min(1,hastings)
  if(alpha>Us[1,l]){
    return(X_new)
  }
  return(this_X_current)
}

shift_update <- function(this_X_current,X_pool,Y,t,l,l_current,F,delta,Us,l_new){
  # shift update: update l then x
  # update of l UAR(1:L) done above since independent of l
  X_new <- this_X_current + F %*% (X_pool[t-1,l_new[l],]-X_pool[t-1,l_current,]) # faster abit
  lnum <- sum(dpois(Y[t,],delta*abs(X_new),log=TRUE))
  ldenom <- sum(dpois(Y[t,],delta*abs(this_X_current),log=TRUE))
  lhastings <- lnum - ldenom
  hastings <- exp(lhastings)
  alpha <- min(1,hastings)
  if(alpha>Us[2,l]){
    return(list(X_new=X_new,l_new=l_new[l]))
  } else {
    return(list(X_new=this_X_current,l_new=l_current))
  }
}

forward_pool <- function(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,delta){
  # metropolis
  # sequential update
  X_pool <- array(0,dim=c(T,L,dim))
  acceptance_rate <- matrix(0,nrow=T,ncol=2)
  ar_accept <- 0
  sh_accept <- 0
  
  # start
  es <- runif(L,0.1,0.4)
  es_2 <- es^2
  zs <- matrix(rnorm(L*dim),ncol=dim,nrow=L)
  Us <- runif(L)
  
  # sample the index for the current states 
  k <- sample(1:L,T,replace=TRUE)
  
  # place the first current state 
  X_pool[1,k[1],] <- X_current[1,]
  
  # reversed transition down from k[1] to 1
  if(k[1]>1){
    for(l in (k[1]-1):1){
      if(l %% 2 == 0){ # if even do usual
        # autoregressive update 
        X_new <- mu_init + sqrt(1-es_2[l])*(X_pool[1,l+1,]-mu_init) + es[l]*sigma_init_L%*%zs[l,]
        # hastings ratio
        lnum <- sum(dpois(Y[1,],delta*abs(X_new),log=TRUE))
        ldenom <- sum(dpois(Y[1,],delta*abs(X_pool[1,l+1,]),log=TRUE))
        lhastings <- lnum - ldenom
        hastings <- exp(lhastings)
        
        alpha <- min(1,hastings)
        if(alpha>Us[l]){
          X_pool[1,l,] <- X_new
          ar_accept <- ar_accept + 1
        } else {
          X_pool[1,l,] <- X_pool[1,l+1,]
        }
      } else { # if odd do flip 
        X_new <- -1*X_pool[1,l+1,]
        X_pool[1,l,] <- X_new
      }

    } 
  }
  
  # forward transition 
  if(k[1]<L){
    for(l in (k[1]+1):L){
      if(l %% 2 == 0){ # if even do flip 
        X_new <- -1*X_pool[1,l-1,]
        X_pool[1,l,] <- X_new
      } else { # if odd do usual
        # autoregressive update 
        X_new <- mu_init + sqrt(1-es_2[l])*(X_pool[1,l-1,]-mu_init) + es[l]*sigma_init_L%*%zs[l,]
        # hastings ratio
        lnum <- sum(dpois(Y[1,],delta*abs(X_new),log=TRUE))
        ldenom <- sum(dpois(Y[1,],delta*abs(X_pool[1,l-1,]),log=TRUE))
        lhastings <- lnum - ldenom
        hastings <- exp(lhastings)
        
        alpha <- min(1,hastings)
        if(alpha>Us[l]){
          X_pool[1,l,] <- X_new
          ar_accept <- ar_accept + 1
        } else {
          X_pool[1,l,] <- X_pool[1,l-1,]
        }
      }
    }
  }
  acceptance_rate[1,] <- c(ar_accept,sh_accept)/(L-1)
  
  # then for t>1
  for(t in 2:T){
    this_X_current <- X_current[t,]
    # stochastic initialisation of l 
    diff <- this_X_current - t(X_pool[t-1,,] %*% F) # for symmetric F
    l_lprob <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    prob <- exp(l_lprob-max(l_lprob)); prob <- prob/sum(prob)
    l_original <- sample(1:L,1,prob=prob)
    
    es <- runif(L,0.1,0.4) # hardcoded
    zs <- matrix(rnorm(L*dim),nrow=L,ncol=dim) 
    Us <- matrix(runif(L*2),nrow=2,ncol=L)
    l_new <- sample(1:L,L,replace=TRUE)
    
    # precompute to save cost
    this_X_pool <- matrix(logical(0),nrow=L,ncol=dim) # for this timestep only 
    this_X_pool[k[t],] <- this_X_current
    mu <- X_pool[t-1,,] %*% F # symmetric F 
    sq_term <- sqrt(1-es^2)
    last <- es*zs%*%sigma_U

    ar_accept <- 0 # to measure acceptance rates 
    sh_accept <- 0 
    
    # reversed transition: shift then AR
    l_current <- l_original
    if(k[t]>1){
      this_X_current <- X_current[t,] # now only used to carry information between AR and shift update
      for(l in (k[t]-1):1){
        if(l %% 2 == 0){ # even then do usual
          # Shift then AR 
          # SHIFT HERE 
          shift_out <- shift_update(this_X_current,X_pool,Y,t,l,l_current,F,delta,Us,l_new)
          # count and then update
          if(all(this_X_current!=shift_out$X_new)){
            sh_accept <- sh_accept + 1
          }
          this_X_current <- shift_out$X_new
          l_current <- shift_out$l_new

          # AUTOREGRESSIVE HERE         
          X_new <- autoregressive_update(this_X_current,mu,Y,t,l,l_current,F,delta,sq_term,Us,last)
          # count and then update
          if(all(this_X_current!=X_new)){
            ar_accept <- ar_accept + 1 
          }
          this_X_current <- X_new 
        } else { # odd do flip
          this_X_current <- -1*this_X_current
          if(l_current %% 2 == 0){ # l is even, then -1
            l_current <- l_current - 1
          } else {
            l_current <- l_current + 1
          }
        }

        # set pool
        this_X_pool[l,] <- this_X_current
      } 
    }
    
    # forward transition: AR then shift 
    l_current <- l_original
    if(k[t]<L){
      this_X_current <- X_current[t,] # now only used to carry information between AR and shift update
      for(l in (k[t]+1):L){
        if(l %% 2 == 0){ # even do flip
          this_X_current <- -1*this_X_current
          if(l_current %% 2 == 0){ # l is even, then -1
            l_current <- l_current - 1
          } else {
            l_current <- l_current + 1
          }
        } else { # odd do usual
          # AR then shift
          # AUTOREGRESSIVE HERE 
          X_new <- autoregressive_update(this_X_current,mu,Y,t,l,l_current,F,delta,sq_term,Us,last)
          # count and update
          if(all(this_X_current!=X_new)){
            ar_accept <- ar_accept + 1 
          }
          this_X_current <- X_new 

          shift_out <- shift_update(this_X_current,X_pool,Y,t,l,l_current,F,delta,Us,l_new)
          # count and then update
          if(all(this_X_current!=shift_out$X_new)){
            sh_accept <- sh_accept + 1
          }
          this_X_current <- shift_out$X_new
          l_current <- shift_out$l_new
        }
        # set pool
        this_X_pool[l,] <- this_X_current
      } 
    }
    
    acceptance_rate[t,] <- c(ar_accept,sh_accept)/(L-1)
    X_pool[t,,] <- this_X_pool
    X_current[t,] <- this_X_current
  } 
  return(list(X_pool=X_pool,acceptance_rate=acceptance_rate))
}

backward_sampling <- function(X_pool,L,T,dim,F,sigma_inv){
  # backward sampling
  X_new <- matrix(0,nrow=T,ncol=dim)
  X_new[T,] <- X_pool[T,sample(1:L,1),] # since UAR on 1 to L @T 
  for(t in (T-1):1){
    # batch calculate the transition probability to the chosen
    diff <- X_new[t+1,] - t(X_pool[t,,] %*% F) # also only for symmetric F 
    lprob <- diag(-1/2*t(diff)%*%sigma_inv%*%diff)
    prob <- exp(lprob-max(lprob)); prob <- prob/sum(prob)
    X_new[t,] <- X_pool[t,sample(1:L,1,prob=prob),]
  }
  return(X_new)
}



lambdaModel2 <- function(ssm,N,L,init=NULL,seed=NULL){
  # sampling epsilon values instead
  # setup
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
  # for backward compatibility 
  if(is.null(ssm$sigma_init_L)){
    print('old ssm object, computing Choleskey for sigma_init')
    sigma_init_L <- t(chol(sigma_init))
  } else {
    sigma_init_L <- ssm$sigma_init_L
  }
  
  if(is.null(ssm$sigma_inv)){
    print('old ssm object, computing inverse for sigma')
    sigma_inv <- chol2inv(sigma_U)
  } else {
    sigma_inv <- ssm$sigma_inv
  }
  
  # X_sample <- array(0,dim=c(N+1,T,dim))
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
  
  acceptance_rate <- array(0,dim=c(2*N,T,2)) # 2 for measuring only the autocorrelation and shift updates
  pb <- txtProgressBar(min=0,max=2*N,title="ehmm",style=3)
  
  for(i in seq(2,2*N,2)){
    # forward sequence
    pool_out <- forward_pool(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,delta)
    
    X_pool <- pool_out$X_pool
    acceptance_rate[i-1,,] <- pool_out$acceptance_rate
    
    X_current <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i,,] <- X_current
    
    # reversed sequence
    pool_out <- forward_pool(X_current[seq(T,1,-1),],Y[seq(T,1,-1),],T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,delta)

    X_pool <- pool_out$X_pool
    acceptance_rate[i,seq(T,1,-1),] <- pool_out$acceptance_rate

    X_current[seq(T,1,-1),] <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i+1,,] <- X_current
    
    setTxtProgressBar(pb, i)
  }
  return(list(X_sample=X_sample[-1,,],N=N,L=L,init=init,seed=seed,X_pool=X_pool,
              acceptance_rate=acceptance_rate*2))
}
