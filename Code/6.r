data <- read.table("../SAT.dat",header = TRUE)
tausq0 <- 2.57
theta.guess <- c(8,8,7,7,7,7,8,8)

gibbs.norm.hier.new = function(y, sigsq, tau0, niter,theta, mu=mean(theta)){
  ###################################################
  # Gibbs sampler for a simplified normal h-model
  ###################################################
	m <- length(theta)
	params <- matrix(NA,niter,m)
	means  <- matrix(NA,niter,1)

  colnames(params) <- paste("theta	",1:m,sep="")
  colnames(means)  <- "mu"

  ## storing initial values
  params[1,] <- theta
  means[1]  <- mu

  for (i in 2:niter){

    # sampling mu's
	# cat(params[i-1,], "\t", tau0, "\t", m, "\n")
	mu <- rnorm(1, mean(params[i-1,]), tau0/m)

	# sampling means
    vstar <- 1/(1/sigsq+1/tau0)
    mstar <- (y/sigsq+mu/tau0)*vstar
    theta <- rnorm(m,mstar,sqrt(vstar))

    ## storing current values
    params[i,] <- theta
    means[i] <- mu
    cat(i,theta,"\n")
  }
 
  return(list(theta=params,mu=means))
}

samp.mutau<-gibbs.norm.hier.new(y=data$esteffect, sigsq=data$seeffect, tausq0, niter=1000,theta=theta.guess)
# warnings()

theta.new.samp<-matrix(NA,nrow=ngrid,ncol=8)
mu.new.samp<-rep(NA,ngrid)
new.samp<-matrix(NA,nrow=ngrid,ncol=9)
start<-c(8,8,7,7,7,7,7,8,8)
#start<-c(8,10,10,10,10,10,10,10,10)  #different starting value
new.sample<-function(tausq0,y,sigma){
  param<-matrix(NA,ngrid,9)
  colnames(param)<-c("mu","theta1","theta2","theta3","theta4","theta5","theta6","theta7","theta8")

  xi<-rep(NA,8)
  gamma<-rep(NA,8)
  param[1,]<-start
  for (i in 2:ngrid){
    param[i,1]<-rnorm(1,sum(param[i-1,2:9])/8,sqrt(tausq0/8))
    for (j in 1:8){
      gamma[j]<-1/(1/tausq0+1/sigma[j]^2)
      xi[j]<-(param[i-1,1]/tausq0+y[j]/sigma[j]^2)*gamma[j]
      param[i,j+1]<-rnorm(1,xi[j],sqrt(gamma[j]))
    }
  }
  return(param)
}
new.samp<-new.sample(tausq0,y,sigma)
t(apply(new.samp,2,quantile,c(0.025,0.5,0.975)))
apply(new.samp,2,mean)

par(mfrow=c(2,4))
for (i in 2:9){
  plot(new.samp[,i],type="l")
}

par(mfrow=c(2,4))
for (i in 2:9){
  acf(new.samp[,i])
}
t(apply(new.samp,2,quantile,c(0.025,0.5,0.975)))
apply(new.samp,2,mean)

## Trying different starting value
start<-c(8,10,10,10,10,10,10,10,10)
