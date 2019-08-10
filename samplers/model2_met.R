metModel2 <- function(ssm,N,es,init=NULL,thin.factor=10,seed=NULL){
  require(MASS)
  # extract ssm values
  delta <- ssm$delta
  F <- ssm$F
  T <- ssm$T
  Q <- ssm$Q
  mu_init <- ssm$mu_init
  sigma_init <- ssm$sigma_init
  # sigma <- ssm$G %*% ssm$Q %*% t(ssm$G)
  sigma <- ssm$sigma
  dim <- ssm$dim
  Y <- ssm$Y

  pre_mu_1 <- solve(F %*% F + solve(sigma_init)%*%sigma) %*% F
  sigma_1 <- solve(F %*% (solve(sigma) %*% F) + solve(sigma_init))
  pre_mu_j <- solve(F %*% F + diag(1,dim)) %*% F
  sigma_j <- solve(F %*% (solve(sigma) %*% F) + solve(sigma))
  pre_mu_T <- F
  sigma_T <- sigma
  
  sigma_1_L <- t(chol(sigma_1))
  sigma_j_L <- t(chol(sigma_j))
  sigma_T_L <- t(chol(sigma_T))
  
  ### begin
  X_sample <- array(0,dim=c(N/thin.factor,T,dim))
  if(is.null(init)){
    if(!is.null(seed)){
      set.seed(seed)
      # print(seed)
    }
    X_sample[1,,] <- mvrnorm(T,mu_init,sigma_init) 
  } else {
    X_sample[1,,] <- init
  }
  X_current <- X_sample[1,,]
  pb <- txtProgressBar(min=0,max=N,title="MH-MCMC",style=3)
  for(i in 2:N){
    # print(paste('progress: ',i/N*100,'%',sep=''))
    setTxtProgressBar(pb, i)
    U <- runif(T)
    for(j in 1:T){
      # find the mu 
      if(j==1) {
        # mu_j <- pre_mu_1 %*% X_sample[i-1,2,]
        mu_j <- pre_mu_1 %*% X_current[2,]
        L <- sigma_1_L
      } else if(j==T) {
        mu_j <- pre_mu_T %*% X_current[T-1,]
        L <- sigma_T_L
      } else {
        # mu_j <- pre_mu_j %*% (X_sample[i,j-1,]+X_sample[i-1,j+1,])
        mu_j <- pre_mu_j %*% (X_current[j-1,]+X_current[j+1,])
        L <- sigma_j_L
      }
      
      # autoregressive update
      z <- rnorm(dim)
      e <- es[i%%length(es)+1] # alternate
      x_j <- mu_j+sqrt(1-e^2)*(X_current[j,]-mu_j)+e*L%*%z
      
      # transition probability
      # lnum <- sum(dpois(Y[j,],exp(c+delta*x_j),log=TRUE))
      # ldenom <- sum(dpois(Y[j,],exp(c+delta*X_current[j,]),log=TRUE))
      lnum <- sum(dpois(Y[j,],delta*abs(x_j),log=TRUE))
      ldenom <- sum(dpois(Y[j,],delta*abs(X_current[j,]),log=TRUE))
      hastings <- exp(lnum-ldenom)
      alpha <- min(1,hastings)
      if(alpha>U[j]){
        X_current[j,] <- x_j
      }
    }
    
    if(i%%thin.factor==0){
      X_sample[i/thin.factor,,] <- X_current
    }
  }
  close(pb)
  return(list(X_sample=X_sample,N=N,thin.factor=thin.factor,es=es,init=init,seed=seed))
}

