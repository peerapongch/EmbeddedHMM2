autoregressive_update <- function(this_X_current,mu,Y,t,l,l_current,F,c,delta,sq_term,Us,last){
  mu_current <- mu[l_current,]
  X_new <- mu_current + sq_term[l]*(this_X_current-mu_current)+last[l,]
  # hastings ratio
  lnum <- sum(dpois(Y[t,],exp(c+delta*X_new),log=TRUE))
  ldenom <- sum(dpois(Y[t,],exp(c+delta*this_X_current),log=TRUE))
  lhastings <- lnum - ldenom
  hastings <- exp(lhastings)
  alpha <- min(1,hastings)
  if(alpha>Us[1,l]){
    return(X_new)
  }
  return(this_X_current)
}

shift_update <- function(this_X_current,X_pool,Y,t,l,l_current,F,c,delta,Us,l_new){
  # shift update: update l then x
  # update of l UAR(1:L) done above since independent of l
  X_new <- this_X_current + F %*% (X_pool[t-1,l_new[l],]-X_pool[t-1,l_current,]) # faster abit
  lnum <- sum(dpois(Y[t,],exp(c+delta*X_new),log=TRUE))
  ldenom <- sum(dpois(Y[t,],exp(c+delta*this_X_current),log=TRUE))
  lhastings <- lnum - ldenom
  hastings <- exp(lhastings)
  alpha <- min(1,hastings)
  if(alpha>Us[2,l]){
    return(list(X_new=X_new,l_new=l_new[l]))
  } else {
    return(list(X_new=this_X_current,l_new=l_current))
  }
}

forward_pool <- function(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,c,delta){
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
  
  # transition down from k[1] to 1
  if(k[1]>1){
    for(l in (k[1]-1):1){
      # autoregressive update 
      X_new <- mu_init + sqrt(1-es_2[l])*(X_pool[1,l+1,]-mu_init) + es[l]*sigma_init_L%*%zs[l,]
      # hastings ratio
      lnum <- sum(dpois(Y[1,],exp(c+delta*X_new),log=TRUE))
      ldenom <- sum(dpois(Y[1,],exp(c+delta*X_pool[1,l+1,]),log=TRUE))
      lhastings <- lnum - ldenom
      hastings <- exp(lhastings)
      
      alpha <- min(1,hastings)
      if(alpha>Us[l]){
        X_pool[1,l,] <- X_new
        ar_accept <- ar_accept + 1
      } else {
        X_pool[1,l,] <- X_pool[1,l+1,]
      }
    } 
  }
  
  if(k[1]<L){
    for(l in (k[1]+1):L){
      # autoregressive update 
      X_new <- mu_init + sqrt(1-es_2[l])*(X_pool[1,l-1,]-mu_init) + es[l]*sigma_init_L%*%zs[l,]
      # hastings ratio
      lnum <- sum(dpois(Y[1,],exp(c+delta*X_new),log=TRUE))
      ldenom <- sum(dpois(Y[1,],exp(c+delta*X_pool[1,l-1,]),log=TRUE))
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
        # Shift then AR 
        # SHIFT HERE 
        shift_out <- shift_update(this_X_current,X_pool,Y,t,l,l_current,F,c,delta,Us,l_new)
        # count and then update
        if(all(this_X_current!=shift_out$X_new)){
          sh_accept <- sh_accept + 1
        }
        this_X_current <- shift_out$X_new
        l_current <- shift_out$l_new

        # AUTOREGRESSIVE HERE         
        X_new <- autoregressive_update(this_X_current,mu,Y,t,l,l_current,F,c,delta,sq_term,Us,last)
        # count and then update
        if(all(this_X_current!=X_new)){
          ar_accept <- ar_accept + 1 
        }
        this_X_current <- X_new 
        
        # set pool
        this_X_pool[l,] <- this_X_current
      } 
    }
    
    # forward transition: AR then shift 
    l_current <- l_original
    if(k[t]<L){
      this_X_current <- X_current[t,] # now only used to carry information between AR and shift update
      for(l in (k[t]+1):L){
        # AR then shift
        # AUTOREGRESSIVE HERE 
        X_new <- autoregressive_update(this_X_current,mu,Y,t,l,l_current,F,c,delta,sq_term,Us,last)
        # count and update
        if(all(this_X_current!=X_new)){
          ar_accept <- ar_accept + 1 
        }
        this_X_current <- X_new 

        shift_out <- shift_update(this_X_current,X_pool,Y,t,l,l_current,F,c,delta,Us,l_new)
        # count and then update
        if(all(this_X_current!=shift_out$X_new)){
          sh_accept <- sh_accept + 1
        }
        this_X_current <- shift_out$X_new
        l_current <- shift_out$l_new

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

makeSigma <- function(rho,dim){
  sigma <- matrix(rho,ncol=dim,nrow=dim)
  diag(sigma) <- 1
  return(sigma)
}

makeSigma_init <- function(rho,phi){
  v <- 1/sqrt((1-phi^2))
  sigma_init <- rho * v %*% t(v)
  diag(sigma_init) <- diag(sigma_init)/rho
  return(sigma_init)
}

param_lprob <- function(param_formed,X_current,Y,T){
  # form
  delta <- param_formed$delta
  c <- param_formed$c
  F <- param_formed$F
  phi <- diag(F)
  sigma_init <- param_formed$sigma_init
  sigma <- param_formed$sigma
  # then lots of precomputation 
  sigma_init_inv <- param_formed$sigma_init_inv
  sigma_inv <- param_formed$sigma_inv

  X_current_t <- t(X_current)
  mu <- phi * X_current_t
  diff <- X_current_t[,2:T] - mu[,1:(T-1)]
  # p(x_1:T|param)
  lprob1 <- -1/2*(log(det(sigma_init)) + (T-1)*log(det(sigma)) + t(X_current[1,])%*%sigma_init_inv%*%X_current[1,] + sum(diag(t(diff)%*%sigma_inv%*%diff)))
  # p(y_1:T|x_1:T,param)
  lambda <- as.vector(exp(c+delta*t(X_current)))
  lprob2 <- sum(dpois(as.vector(t(Y)),lambda,log=TRUE))
  
  return(lprob1+lprob2)
}

update_shsc_param <- function(param_current,param_formed_current,param_lprob_current,X_current,Y,shsc_scale,dim,T){
  # propose
  param_new <- param_current
  z <- rnorm(2)*shsc_scale
  param_new[3:4] <- param_new[3:4] + z

  check_domain_delta <- (param_new[3] > 0)

  if(check_domain_delta){
    # form 
    param_formed_new <- param_formed_current
    param_formed_new$delta <- rep(param_new[3],dim)
    param_formed_new$c <- rep(param_new[4],dim)
    X_new <- (X_current*param_formed_current$delta[1]+param_formed_current$c[1]-param_formed_new$c[1])/param_formed_new$delta[1]
    
    # compute prob
    param_lprob_new <- param_lprob(param_formed_new,X_new,Y,T)
    N <- T*dim
    log_acceptance_prob <- param_lprob_new-param_lprob_current+N*(log(param_formed_current$delta[1])-log(param_formed_new$delta[1]))

    if(log(runif(1))<log_acceptance_prob){
      # print('ACCEPTED!!')
      param_current <- param_new
      # print(param_current)
      X_current <- X_new
      param_lprob_current <- param_lprob_new
      param_formed_current <- param_formed_new
    }
  } # else automatic out and do nothing but return the current values

  return(list(param_new=param_current,param_lprob_new=param_lprob_current,param_formed_new=param_formed_current,X_new=X_current))
}

update_rw_param <- function(param_current,param_formed_current,param_lprob_current,X_current,Y,rw_scale,dim){
  # propose
  param_new <- param_current + rw_scale*rnorm(4)

  check_domain_phi <- (param_new[1] > 0) && (param_new[1] < 1)
  check_domain_rho <- (param_new[2] > 0) && (param_new[2] < 1)
  check_domain_delta <- (param_new[3] > 0)

  if(check_domain_phi && check_domain_rho && check_domain_delta){
    # form 
    F <- diag(rep(param_new[1],dim))
    delta_new <- rep(param_new[3],dim)
    c_new <- rep(param_new[4],dim)
    sigma_init <- makeSigma_init(rep(param_new[2],dim),rep(param_new[1],dim))
    sigma <- makeSigma(rep(param_new[2],dim),dim)

    # then lots of precomputation 
    sigma_init_U <- chol(sigma_init)
    sigma_init_L <- t(sigma_init_U)
    sigma_init_inv <- chol2inv(sigma_init_U)
    sigma_U <- chol(sigma)
    sigma_L <- t(sigma_U)
    sigma_inv <- chol2inv(sigma_U)
    param_formed_new <- list(F=F,sigma_init=sigma_init,sigma=sigma,sigma_init_inv=sigma_init_inv,sigma_inv=sigma_inv,
      sigma_init_U = sigma_init_U, sigma_init_L = sigma_init_L, sigma_U = sigma_U, sigma_L = sigma_L,
      delta=delta_new,c=c_new)
    
    # compute prob
    param_lprob_new <- param_lprob(param_formed_new,X_current,Y,T)

    if(log(runif(1))<param_lprob_new-param_lprob_current){
      # print('ACCEPTED!!')
      param_current <- param_new
      # print(param_current)
      param_lprob_current <- param_lprob_new
      param_formed_current <- param_formed_new
    }
  } # else automatic out and do nothing but return the current values

  return(list(param_new=param_current,param_lprob_new=param_lprob_current,param_formed_new=param_formed_current))
}

lambdaModel1_param <- function(ssm,N,L,init=NULL,seed=NULL,N.mcmc.param=20,rw_scale=c(0.1,0.1,0.1,0.1),shsc_scale=c(0.01,0.1),rho.init=NULL,phi.init=NULL,delta.init=NULL,c.init=NULL,checkpoint.name=NULL){
  # sampling epsilon values instead
  # setup
  require(MASS)
  if(!is.null(seed)){
    set.seed(seed)
  }
  #poisson observation and gaussian latent process 
  Y <- ssm$Y
  dim <- ssm$dim
  T <- ssm$T
  mu_init <- ssm$mu_init # later include in the sampling step 

  if(is.null(phi.init)){
    phi.init <- runif(1)
  }

  if(is.null(rho.init)){
    rho.init <- runif(1)
  }

  if(is.null(delta.init)){
    delta.init <- runif(1)
  }

  if(is.null(c.init)){
    c.init <- rnorm(1)
  }  

  # parameters to be sampled
  param_sample <- matrix(logical(0),nrow=N.mcmc.param*2*N+1,ncol=4) # order is phi, rho; delta
  param_current <- c(phi.init,rho.init,delta.init,c.init)
  print(param_current)
  param_sample[1,] <- param_current

  # form
  delta <- rep(delta.init,dim)
  c <- rep(c.init,dim)
  F <- diag(rep(phi.init,dim))
  sigma_init <- makeSigma_init(rep(rho.init,dim),rep(phi.init,dim))
  sigma <- makeSigma(rep(rho.init,dim),dim)
  # then lots of precomputation 
  sigma_init_U <- chol(sigma_init)
  sigma_init_L <- t(sigma_init_U)
  sigma_init_inv <- chol2inv(sigma_init_U)
  sigma_U <- chol(sigma)
  sigma_L <- t(sigma_U)
  sigma_inv <- chol2inv(sigma_U)

  param_formed_current <- list(F=F,sigma_init=sigma_init,sigma=sigma,sigma_init_inv=sigma_init_inv,sigma_inv=sigma_inv,
    sigma_init_U = sigma_init_U, sigma_init_L = sigma_init_L, sigma_U = sigma_U, sigma_L = sigma_L,
    delta=delta,c=c)
  
  # X_sample <- array(0,dim=c(N+1,T,dim))
  X_sample <- array(0,dim=c(2*N+1,T,dim))
  if(is.null(init)){
    X_sample[1,,] <- mvrnorm(T,mu_init,sigma_init)
  } else {
    X_sample[1,,] <- init 
  }
  X_current <- X_sample[1,,]
  
  acceptance_rate <- array(0,dim=c(2*N,T,2)) # 2 for measuring only the autocorrelation and shift updates
  pb <- txtProgressBar(min=0,max=2*N,title="ehmm",style=3)
  rw_acceptance_rate <- 0
  shsc_acceptance_rate <- 0

  for(i in seq(2,2*N,2)){
    if(!is.null(checkpoint.name)){
      if((i == N/2) || (i == N) || (i == N/4*3)){
        save(X_sample,param_sample,file=checkpoint.name)
      }
    }
    # forward sequence
    pool_out <- forward_pool(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,c,delta)
    
    X_pool <- pool_out$X_pool
    acceptance_rate[i-1,,] <- pool_out$acceptance_rate
    
    X_current <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i,,] <- X_current

    # parameter sampling 1
    # sample
    param_lprob_current <- param_lprob(param_formed_current,X_current,Y,T)
    start_index <- N.mcmc.param*(i-1) - (N.mcmc.param-1)
    for(j in seq(1,N.mcmc.param,2)){
      # all together 
      param_update_out <- update_rw_param(param_current,param_formed_current,param_lprob_current,X_current,Y,rw_scale,dim)
      if(any(param_current != param_update_out$param_new)){
        rw_acceptance_rate <- rw_acceptance_rate + 1
      }
      param_formed_current <- param_update_out$param_formed_new
      param_current <- param_update_out$param_new
      param_lprob_current <- param_update_out$param_lprob_new
      param_sample[start_index+j,] <- param_current

      # then shsc update 
      param_update_out <- update_shsc_param(param_current,param_formed_current,param_lprob_current,X_current,Y,shsc_scale,dim,T)
      if(any(param_current != param_update_out$param_new)){
        shsc_acceptance_rate <- shsc_acceptance_rate + 1
      }
      param_formed_current <- param_update_out$param_formed_new
      param_current <- param_update_out$param_new
      X_current <- param_update_out$X_new
      param_lprob_current <- param_update_out$param_lprob_new
      param_sample[start_index+j+1,] <- param_current
    }

    # update parameters
    delta <- param_formed_current$delta
    c <- param_formed_current$c
    F <- param_formed_current$F
    sigma_init <- param_formed_current$sigma_init
    sigma <- param_formed_current$sigma
    sigma_init_U <- param_formed_current$sigma_init_U
    sigma_init_L <- param_formed_current$sigma_init_L
    sigma_init_inv <- param_formed_current$sigma_init_inv
    sigma_U <- param_formed_current$sigma_U
    sigma_L <- param_formed_current$sigma_L
    sigma_inv <- param_formed_current$sigma_inv
    
    # reversed sequence
    pool_out <- forward_pool(X_current[seq(T,1,-1),],Y[seq(T,1,-1),],T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,c,delta)

    X_pool <- pool_out$X_pool
    acceptance_rate[i,seq(T,1,-1),] <- pool_out$acceptance_rate

    X_current[seq(T,1,-1),] <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i+1,,] <- X_current

    # parameter sampling 2 
    # sample
    param_lprob_current <- param_lprob(param_formed_current,X_current,Y,T)
    start_index <- N.mcmc.param*(i) - (N.mcmc.param-1)
    for(j in seq(1,N.mcmc.param,2)){
      # all together 
      param_update_out <- update_rw_param(param_current,param_formed_current,param_lprob_current,X_current,Y,rw_scale,dim)
      if(any(param_current != param_update_out$param_new)){
        rw_acceptance_rate <- rw_acceptance_rate + 1
      }
      param_formed_current <- param_update_out$param_formed_new
      param_current <- param_update_out$param_new
      param_lprob_current <- param_update_out$param_lprob_new
      param_sample[start_index+j,] <- param_current

      # then shsc update 
      param_update_out <- update_shsc_param(param_current,param_formed_current,param_lprob_current,X_current,Y,shsc_scale,dim,T)
      if(any(param_current != param_update_out$param_new)){
        shsc_acceptance_rate <- shsc_acceptance_rate + 1
      }
      param_formed_current <- param_update_out$param_formed_new
      param_current <- param_update_out$param_new
      X_current <- param_update_out$X_new
      param_lprob_current <- param_update_out$param_lprob_new
      param_sample[start_index+j+1,] <- param_current
    }

    # update parameters
    delta <- param_formed_current$delta
    c <- param_formed_current$c
    F <- param_formed_current$F
    sigma_init <- param_formed_current$sigma_init
    sigma <- param_formed_current$sigma
    sigma_init_U <- param_formed_current$sigma_init_U
    sigma_init_L <- param_formed_current$sigma_init_L
    sigma_init_inv <- param_formed_current$sigma_init_inv
    sigma_U <- param_formed_current$sigma_U
    sigma_L <- param_formed_current$sigma_L
    sigma_inv <- param_formed_current$sigma_inv
    
    setTxtProgressBar(pb, i)
  }
  rw_acceptance_rate <- rw_acceptance_rate/(N.mcmc.param*N)
  shsc_acceptance_rate <- shsc_acceptance_rate/(N.mcmc.param*N)

  return(list(X_sample=X_sample[-1,,],param_sample=param_sample,N=N,L=L,init=init,seed=seed,X_pool=X_pool, acceptance_rate=acceptance_rate,rw_acceptance_rate=rw_acceptance_rate,shsc_acceptance_rate=shsc_acceptance_rate))
}
