
gibbs.norm.mixture.new = function(y,n,niter,theta){
  ###################################################
  # Gibbs sampler for mixture of two normals 
  ###################################################
  params <- matrix(NA,niter,5)
  indics <- matrix(NA,niter,n)

  colnames(params) <- c("alpha","mu0","mu1","sigsq0","sigsq1")
  colnames(indics) <- paste("z",1:n,sep="")
  
  ## storing initial values
  params[1,] <- theta
  indics[1,] <- rep(0,n)

  ## gibbs sampler
  for (i in 2:niter){

    ## sampling indicator variables
    p0 <- dnorm(y,theta[2],sqrt(theta[4]))*(1-theta[1])
    p1 <- dnorm(y,theta[3],sqrt(theta[5]))*(theta[1])
    z <- rbinom(n,1,p1/(p1+p0))

    ## calculating statistics from indicator variables
    n0 <- sum(z==0)
    n1 <- sum(z==1)

    ## sampling alpha
    theta[1] <- rbeta(1,n1+1,n0+1)

    ## sampling means
    theta[2] <- rnorm(1,mean(y[z==0]),sqrt(theta[4]/n0))
    theta[3] <- rnorm(1,mean(y[z==1]),sqrt(theta[5]/n1))

    ## sampling variances
    theta[4] <- 1/rgamma(1,n0/2+1,sum((y[z==0]-theta[2])^2)/2)
    theta[5] <- 1/rgamma(1,n1/2+1,sum((y[z==1]-theta[3])^2)/2)

    ## storing current values
    params[i,] <- theta
    indics[i,] <- z
    cat(i,theta,"\n")
  }

  return(list(params=params,indics=indics))
}

gibbs.norm.hier.new = function(y,ybar,m,n,ntot,niter,theta){
  ###################################################
  # Gibbs sampler for a normal hierarchical model 
  ###################################################
  params <- matrix(NA,niter,4)
  means  <- matrix(NA,niter,m)

  colnames(params) <- c("mu0","tausq","sigsq","logpost")
  colnames(means)  <- paste("mu",1:m,sep="")

  ## storing initial values
  params[1,] <- theta
  means[1,]  <- mu

  for (i in 2:niter){

    # sampling mu's
    vstar <- 1/(n/sigsq+1/tausq)
    mstar <- (ybar*n/sigsq+mu0/tausq)*vstar
    mu <- rnorm(m,mstar,sqrt(vstar))

    # sampling mu0
    theta[1] <- rnorm(1,mean(mu),sqrt(theta[2]/m))

    # sampling tausq
    theta[2] <- 1/rgamma(1,(m-1)/2,sum((mu-theta[1])^2)/2)

    # sampling sigsq
    sumsq.y <- 0
    for (j in 1:m) {sumsq.y <- sumsq.y + sum((y[,j]-mu[j])^2)}
    theta[3] <- 1/rgamma(1,ntot/2,sumsq.y/2)

    ## storing current values
    params[i,] <- theta
    means[i,] <- mu
    cat(i,theta,"\n")
  }
 
  return(list(theta=params,mu=means))
}




logpost.norm.hier <- function(ybar,n,mu,theta){
  #############################################################
  # Evaluate the joint posterior for normal hierarchical model
  #############################################################
  return(sum(dnorm(mu,theta[1],sqrt(theta[2]),log=T))+sum(dnorm(ybar,mu,sqrt(theta[3]/n),log=T)))
}

logpost.norm.mixture <- function(Y,Z,theta){
  ########################################################
  # Evaluate the joint posterior for normal mixture model
  ########################################################
  mean <- ifelse(Z==0,theta[2],theta[3])
  sdev <- ifelse(Z==0,sqrt(theta[4]),sqrt(theta[5]))
  return(sum(dbinom(Z,1,theta[1],log=T))+sum(dnorm(Y,mean,sdev,log=T)))
}

gibbs.norm.hier = function(y,ybar,m,n,ntot,niter,start,plot=F,eval.post=F){
  ###################################################
  # Gibbs sampler for a normal hierarchical model 
  ###################################################
  means <- matrix(NA,niter,m)
  params <- matrix(NA,niter,4)

  colnames(means) = paste("mu",1:m,sep="")
  colnames(params) = c("mu0","tausq","sigsq","logpost")

  mu0   <- start[1]
  tausq <- start[2]
  sigsq <- start[3]
  mu    <- rep(0,m)

  logpost <- 0
  if (eval.post) logpost <- logpost.norm.hier(ybar,n,mu,start)

  params[1,] <- c(start,logpost)
  means[1,] <- mu

  #if (plot) {graphics.off(); dev.new(); dev.new()}

  for (i in 2:niter){

    # sampling mu's
    vstar <- 1/(n/sigsq+1/tausq)
    mstar <- (ybar*n/sigsq+mu0/tausq)*vstar
    mu <- rnorm(m,mstar,sqrt(vstar))

    # sampling mu0
    mu0 <- rnorm(1,mean(mu),sqrt(tausq/m))

    # sampling tausq
    tausq <- 1/rgamma(1,(m-1)/2,sum((mu-mu0)^2)/2)

    # sampling sigsq
    sumsq.y <- 0
    for (j in 1:m) {sumsq.y <- sumsq.y + sum((y[,j]-mu[j])^2)}
    sigsq <- 1/rgamma(1,ntot/2,sumsq.y/2)

    ## evaluate the log posterior   
    theta <- c(mu0,tausq,sigsq)
    if (eval.post) logpost <- logpost.norm.hier(ybar,n,mu,theta)

    ## storing current values
    means[i,] <- mu
    params[i,] <- c(theta,logpost)
    cat(i,mu0,tausq,sigsq,logpost,"\n")

    ## plot iterations 
    good <- c(1,2,9)
    plot.time <- plot & ((i<=500 & i%%100==0) | (i>500 & i%%(niter/4)==0))
    if (plot.time){
      if (i<=500) {dev.set(2); plot.ts(params);           dev.set(3); plot.ts(means[,good])} 
      else        {dev.set(2); plot.ts(params[-(1:10),]); dev.set(3); plot.ts(means[-(1:10),good])}
    }
  }
  if (plot) {dev.set(2); boxplot(y)}
  if (plot) {dev.set(3); boxplot(means[-(1:10),]~col(means[-(1:10),]),ylim=range(y))}
 
  return(list(mu=means,theta=params))
}


gibbs.norm.mixture = function(Y,n,niter,start,plot=F,eval.post=F){
  ###################################################
  # Gibbs sampler for mixture of two normals 
  ###################################################
  params <- matrix(NA,niter,6)
  indics <- matrix(NA,niter,n)

  colnames(indics) <- paste("Z",1:n,sep="")
  colnames(params) <- c("alpha","mu0","mu1","sigsq0","sigsq1","logpost")
  

  if (plot) {graphics.off(); dev.new(); dev.new()}

  #initializing parameters and mixture indicators
  alpha  <- start[1]
  mu0    <- start[2]   
  mu1    <- start[3]
  sigsq0 <- start[4]   
  sigsq1 <- start[5]
  Z      <- rep(0,n)

  logpost <- 0
  if (eval.post) logpost <- logpost.norm.mixture(Y,Z,start)

  params[1,] <- c(start,logpost)
  indics[1,] <- Z

  ## gibbs sampler
  for (i in 2:niter){

    ## sampling indicator variables
    p1 <- dnorm(Y,mu1,sqrt(sigsq1))*alpha    
    p0 <- dnorm(Y,mu0,sqrt(sigsq0))*(1-alpha)
    Z <- rbinom(n,1,p1/(p1+p0))

    ## calculating statistics from indicator variables
    n0 <- sum(Z==0)
    n1 <- sum(Z==1)
    mean0 <- mean(Y[Z==0])
    mean1 <- mean(Y[Z==1])

    ## sampling alphas
    alpha <- rbeta(1,n1+1,n0+1)

    ## sampling means
    mu0 <- rnorm(1,mean0,sqrt(sigsq0/n0))
    mu1 <- rnorm(1,mean1,sqrt(sigsq1/n1))

    ## calculating statistics from new means 
    ss0 <- sum((Y[Z==0]-mu0)^2)
    ss1 <- sum((Y[Z==1]-mu1)^2)

    ## sampling variances
    sigsq0 <- 1/rgamma(1,n0/2+1,ss0/2)
    sigsq1 <- 1/rgamma(1,n1/2+1,ss1/2)

    ## evaluate the log posterior   
    theta <- c(alpha,mu0,mu1,sigsq0,sigsq1)
    if (eval.post) logpost <- logpost.norm.mixture(Y,Z,theta)

    ## storing current values
    params[i,] <- c(theta,logpost)
    indics[i,] <- Z
    cat(i,alpha,mu0,mu1,sigsq0,sigsq1,logpost,"\n")

    ## plot iterations 
    good <- seq(1,n,len=3)
    plot.time <- plot & ((i<=500 & i%%100==0) | (i>500 & i%%(niter/4)==0))
    if (plot.time){
      if (i<=500) {dev.set(2); plot.ts(params[1:i,]);   dev.set(3); plot.ts(indics[1:i,good])} 
      else        {dev.set(2); plot.ts(params[101:i,]); dev.set(3); plot.ts(indics[101:i,good])}
    }
  }

  #if (plot) dev.set(3); pairs(params[-(1:100),1:5])

  return(list(params=params,indics=indics))
}
