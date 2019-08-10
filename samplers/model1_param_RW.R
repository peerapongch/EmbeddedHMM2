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

param_rw_lprob <- function(param_formed,X_current,Y,T){
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
  # print(lprob1)
  # p(y_1:T|x_1:T,param)

  lambda <- as.vector(exp(c+delta*t(X_current)))
  # print(delta)
  # print(lambda)
  lprob2 <- sum(dpois(as.vector(t(Y)),lambda,log=TRUE))
  # print(dpois(as.vector(t(Y)),lambda,log=TRUE))

  return(lprob1+lprob2)
}

param_update_rw <- function(param_current,param_formed_current,param_lprob_current,X_current,Y,rw_scale,dim){
  # propose
  param_new <- param_current + rw_scale*rnorm(4)

  check_domain_phi <- (param_new[1] > 0) && (param_new[1] < 1)
  check_domain_rho <- (param_new[2] > 0) && (param_new[2] < 1)
  check_domain_delta <- (param_new[3] > 0)
  check_domain_c <- (param_new[4] > -10) && (param_new[4] < 10)

  if(check_domain_phi && check_domain_rho && check_domain_delta && check_domain_c){
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
    param_lprob_new <- param_rw_lprob(param_formed_new,X_current,Y,T)

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

paramModel1 <- function(ssm,N,X=NULL,seed=NULL,rw_scale=c(0.1,0.1,0.1,0.1),rho.init=NULL,phi.init=NULL,delta.init=NULL,c.init=NULL,checkpoint.name=NULL){
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
  param_sample <- matrix(logical(0),nrow=N+1,ncol=4) # order is phi, rho; delta
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

  param_formed_current <- list(F=F,sigma_init=sigma_init,sigma=sigma,sigma_init_inv=sigma_init_inv,sigma_inv=sigma_inv, sigma_init_U = sigma_init_U, sigma_init_L = sigma_init_L, sigma_U = sigma_U, sigma_L = sigma_L, delta=delta,c=c)
  
  if(is.null(X)){
    X_current <- ssm$X
  } else {
    X_current <- X
  }
  
  pb <- txtProgressBar(min=0,max=N+1,style=3)
  param_acceptance_rate <- 0

  for(i in 1:(N+1)){
    if(!is.null(checkpoint.name)){
      if(i == N){
        save(X_sample,param_sample,file=checkpoint.name)
      }
    }
    param_lprob_current <- param_rw_lprob(param_formed_current,X_current,Y,T)
    param_update_out <- param_update_rw(param_current,param_formed_current,param_lprob_current,X_current,Y,rw_scale,dim)
    if(any(param_current != param_update_out$param_new)){
        param_acceptance_rate <- param_acceptance_rate + 1
      }
    param_formed_current <- param_update_out$param_formed_new
    param_current <- param_update_out$param_new
    param_lprob_current <- param_update_out$param_lprob_new
    param_sample[i,] <- param_current

    setTxtProgressBar(pb, i)
  }
  param_acceptance_rate <- param_acceptance_rate/N

  return(list(param_sample=param_sample,N=N,L=L,seed=seed,param_acceptance_rate=param_acceptance_rate))
}
