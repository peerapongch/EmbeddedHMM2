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

autoregressive_update <- function(this_X_current,mu,Y,t,l,l_current,F,delta,sq_term,Us,last){
# autoregressive_update <- function(this_X_current,X_pool,Y,t,l,l_current,F,sigma_L,c,delta,es,es_2,zs,Us){
  # mu_current <- F %*% X_pool[t-1,l_current,]
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

mcmc_step<- function(X_current,Y,N,es,T,dim,mus,sigmas_L,delta,thin.factor=10,seed=NULL){
  require(MASS)
  
  ##### MHMCMC autoregressive ######
  X_sample <- array(0,dim=c(N/thin.factor,T,dim))
  X_sample[1,,] <- X_current
  if(!is.null(seed)){
    set.seed(seed)
  }
  for(i in 2:N){
    U <- runif(T)
    for(j in 1:T){
      # find the mu 
      if(j==1) {
        # mu_j <- pre_mu_1 %*% X_sample[i-1,2,]
        mu_j <- mus[[1]] %*% X_current[2,]
        L <- sigmas_L[[1]]
      } else if(j==T) {
        mu_j <- mus[[3]] %*% X_current[T-1,]
        L <- sigmas_L[[3]]
      } else {
        # mu_j <- pre_mu_j %*% (X_sample[i,j-1,]+X_sample[i-1,j+1,])
        mu_j <- mus[[2]] %*% (X_current[j-1,]+X_current[j+1,])
        L <- sigmas_L[[2]]
      }
      
      # autoregressive update
      z <- rnorm(dim)
      e <- es[i%%length(es)+1] # alternate
      x_j <- mu_j+sqrt(1-e^2)*(X_current[j,]-mu_j)+e*L%*%z
      
      lnum <- sum(dpois(Y[j,],delta*abs(x_j),log=TRUE))
      ldenom <- sum(dpois(Y[j,],delta*abs(X_current[j,]),log=TRUE))
      hastings <-exp(lnum - ldenom)
      
      alpha <- min(1,hastings)
      if(alpha>U[j]){
        X_current[j,] <- x_j
      } #else unchanged
    }
    
    if(i%%thin.factor==0){
      X_sample[i/thin.factor,,] <- X_current
    }
  }
  return(X_current) # return only the last sample
  # return(list(X_sample=X_sample,N=N,thin.factor=thin.factor,es=es,init=init))
}

combinedModel2 <- function(ssm,N,L,L_particles,es,N.mcmc=10,init=NULL,seed=NULL){
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

  # set up for mcmc
  pre_mu_1 <- solve(F %*% F + solve(sigma_init)%*%sigma) %*% F
  sigma_1 <- solve(F %*% (solve(sigma) %*% F) + solve(sigma_init))
  pre_mu_j <- solve(F %*% F + diag(1,dim(F)[1])) %*% F
  sigma_j <- solve(F %*% (solve(sigma) %*% F) + solve(sigma))
  pre_mu_T <- F
  sigma_T <- sigma
  
  sigma_1_L <- t(chol(sigma_1))
  sigma_j_L <- t(chol(sigma_j))
  sigma_T_L <- t(chol(sigma_T))
  
  mus <- list(pre_mu_1,pre_mu_j,pre_mu_T)
  sigmas_L <- list(sigma_1_L,sigma_j_L,sigma_T_L)

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
  X_sample <- array(0,dim=c(6*N+1,T,dim))
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
  
  acceptance_rate <- array(0,dim=c(6*N,T,2)) # 2 for measuring only the autocorrelation and shift updates
  pb <- txtProgressBar(min=0,max=6*N,title="ehmm",style=3)
  
  for(i in seq(2,6*N,6)){
    # Alternate between ehmm and pgmet 
    # ehmm: forward sequence
    pool_out <- forward_pool(X_current,Y,T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,delta)
    
    X_pool <- pool_out$X_pool
    acceptance_rate[i-1,,] <- pool_out$acceptance_rate
    
    X_current <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i,,] <- X_current
    
    # ehmm: reversed sequence
    pool_out <- forward_pool(X_current[seq(T,1,-1),],Y[seq(T,1,-1),],T,L,dim,mu_init,sigma_init_L,sigma_L,sigma_U,F,sigma_inv,delta)

    X_pool <- pool_out$X_pool
    acceptance_rate[i,seq(T,1,-1),] <- pool_out$acceptance_rate

    X_current[seq(T,1,-1),] <- backward_sampling(X_pool,L,T,dim,F,sigma_inv)
    X_sample[i+1,,] <- X_current

    # pgmet: step 1
    forward_results <- forward_step(X_current,Y,T,L_particles,dim,F,sigma_U,mu_init,sigma_init,delta)
    X_current <- backward_step(forward_results,T,L_particles,F,sigma_inv)
    X_sample[i+2,,] <- X_current

    # pgmet: step 2
    X_current <- mcmc_step(X_current,Y,N.mcmc,es,T,dim,mus,sigmas_L,delta,thin.factor=1)
    X_sample[i+3,,] <- X_current

    # pgmet: step 3 
    forward_results <- forward_step(X_current[seq(T,1,-1),],Y[seq(T,1,-1),],T,L_particles,dim,F,sigma_U,mu_init,sigma_init,delta)
    X_current[seq(T,1,-1),] <- backward_step(forward_results,T,L_particles,F,sigma_inv) # reverse the index when setting
    X_sample[i+4,,] <- X_current 

    # pgmet: step 4 
    X_current <- mcmc_step(X_current,Y,N.mcmc,es,T,dim,mus,sigmas_L,delta,thin.factor=1)
    X_sample[i+5,,] <- X_current

    setTxtProgressBar(pb, i)
  }
  return(list(X_sample=X_sample[-1,,],N=N,L=L,L_particles=L_particles,es=es,init=init,seed=seed,X_pool=X_pool,acceptance_rate=acceptance_rate))
}
