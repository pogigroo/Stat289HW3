data <- read.table("../SAT.dat",header = TRUE)
tausq0 <- 2.57
theta.guess <- c(9,9,9,9,9,9,9,9)

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
		mu <- rnorm(1, mean(params[i-1,]), tau0/m)
		
		# sampling means
		vstar <- 1/(1/sigsq+1/tau0)
		mstar <- (y/sigsq+mu/tau0)*vstar
		theta <- rnorm(m,mstar,sqrt(vstar))
		
		## storing current values
		params[i,] <- theta
		means[i] <- mu
		# cat(i,theta,"\n")
	}
 
  return(list(theta=params,mu=means))
}

### I had to generate alot of samples to get ~1000 good draws
gibbs<-gibbs.norm.hier.new(y=data$esteffect, sigsq=data$seeffect, tausq0, niter=10000,theta=theta.guess)
# warnings()

theta <- gibbs$theta
mu <- gibbs$mu


names.theta <- expression(theta[1],theta[2],theta[3],theta[4],theta[5],theta[6],theta[7],theta[8])


### checking convergence 
par(mfrow=c(5,2), mar=c(3.5,3.5,2,1), mgp=c(2,.65,0), las=1, cex.main=1.5)

plot(mu, type='l', main=expression(mu),ylab = "")
for (i in 1:8) {
	plot(theta[,i],type='l', main=names.theta[i], ylab = "")
}


### throwing away first 100 as burn-in
good <- 101:10000

### autocorrelation of draws in chain 1
par(mfrow=c(5,2), mar=c(3.5,3.5,2,1), mgp=c(2,.65,0), las=1, cex.main=1.5)
for (i in 1:8) {acf(theta[good,i],main="", xlab=paste('Lag ',names.theta[i]),ylab = "")}
acf(mu[good],main="",xlab = paste('Lag ', expression(mu)),ylab = "")


### thinning chains (every 10th iteration)
par(mfrow=c(5,2), mar=c(3.5,3.5,2,1), mgp=c(2,.65,0), las=1, cex.main=1.5)
good <- seq(101,10000,by=10)
theta <- theta[good,] #rbind,theta2[good,])
mu <- mu[good] #rbind,mu2[good,])
for (i in 1:8) {acf(theta[,i] , main="", xlab=names.theta[i] ,ylab = "")}
acf(mu,main="",xlab = expression(mu),ylab = "")

### generate quantiles table
mean.quan <- function(x) c(mean(x),quantile(x, c(0.025,0.5,0.975)))
cbind(apply(theta, 2, mean.quan) , mean.quan(mu))

### Problem 7
sum(theta[1,] > theta[2:8,])/nrow(theta)
