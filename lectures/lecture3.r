############################################################
############## R CODE FOR STAT 289 LECTURE 3 ###############
############################################################

##### reading in Poisson data ####

data <- read.table("planes.txt",header=T)
attach(data)

hist(fatal)



sumfatal <- sum(fatal)
n <- length(fatal)

##### looking at different Gamma priors #####
theta <- ppoints(1000)*30
gammaprior1 <- dgamma(theta,shape=10,rate=10)
gammaprior2 <- dgamma(theta,shape=50,rate=10)
gammaprior3 <- dgamma(theta,shape=100,rate=10)
ylim = range(gammaprior1,gammaprior2,gammaprior3)

plot(theta,gammaprior1,type="l",ylim=ylim,lwd=2)
lines(theta,gammaprior2,col=2,lwd=2)
lines(theta,gammaprior3,col=4,lwd=2)
legend(15,1,c("Gamma(10,10)","Gamma(50,10)","Gamma(100,10)"),col=c(1,2,4),lwd=2)





##### looking at different "non-informative" Gamma priors #####
theta <- ppoints(1000)*30
gammaprior4 <- dgamma(theta,shape=1,rate=1)
gammaprior5 <- dgamma(theta,shape=0.5,rate=0.0001)
gammaprior6 <- dgamma(theta,shape=0.0001,rate=0.0001)
ylim = range(gammaprior4,gammaprior5,gammaprior6)

plot(theta,gammaprior4,type="l",ylim=ylim,lwd=2)
lines(theta,gammaprior5,col=2,lwd=2)
lines(theta,gammaprior6,col=4,lwd=2)
legend(15,1,c("Gamma(1,1)","Gamma(0.5,0)","Gamma(0,0)"),col=c(1,2,4),lwd=2)



##### looking at posterior for all these different priors #####
theta <- ppoints(1000)*50
gammaposterior1 <- dgamma(theta,shape=(sumfatal+10),rate=(n+10))
gammaposterior2 <- dgamma(theta,shape=(sumfatal+50),rate=(n+10))
gammaposterior3 <- dgamma(theta,shape=(sumfatal+100),rate=(n+10))
gammaposterior4 <- dgamma(theta,shape=(sumfatal+1),rate=(n+1))
gammaposterior5 <- dgamma(theta,shape=(sumfatal+0.5),rate=(n+0.0001))
gammaposterior6 <- dgamma(theta,shape=(sumfatal+0.0001),rate=(n+0.0001))

ylim = range(gammaposterior1,gammaposterior2,gammaposterior3,
             gammaposterior4,gammaposterior5,gammaposterior6)

par(mfrow=c(2,1),mar=c(3,3,2,1),las=1)
plot(theta,gammaposterior1,type="l",ylim=ylim,lwd=2)
lines(theta,gammaposterior2,col=2,lwd=2)
lines(theta,gammaposterior3,col=3,lwd=2)
lines(theta,gammaposterior4,col=4,lwd=2)
lines(theta,gammaposterior5,col=5,lwd=2)
lines(theta,gammaposterior6,col=6,lwd=2)
abline(v=mean(fatal),col="brown",lwd=2)

hist(fatal,xlim=c(0,50))#,ylim=c(minplot,maxplot),prob=T)
abline(v=mean(fatal),col="brown",lwd=2)



###############################################################
## Normal Simulations

norm.post = function(mu0,tausq0,ybar,n,sigmasq){
  ################################################
  ## Function to summarize normal posteriors
  ################################################
  postvar = 1/(1/tausq0+n/sigmasq)
  postmean = postvar*(mu0/tausq0+ybar*n/sigmasq)
  cat("\n",format(round(mu0,3),nsmall=3),format(round(tausq0,3),nsmall=3),
      "| ",format(round(ybar,3),nsmall=3),format(round(sigmasq/n,3),nsmall=3),
      "| ",format(round(postmean,3),nsmall=3),format(round(postvar,3),nsmall=3),"\n")
  return(list(postmean=postmean,postvar=postvar))
}


##### simulating some Normal data ####
mu.true <- 3
sigsq <- 5
n <- 100
data <- rnorm(n,mean=mu.true,sd=sqrt(sigsq))
par(mfrow=c(2,1),mar=c(3,3,2,1),las=1)
hist(data,main=paste("Normal Data (n=",length(data), ", ybar=",round(mean(data),3),")",sep=""))



## posterior distribution for priors with mu0 = 1 and different variances
mu0 <- 1
tausqs <- c(0.1,1,10)

mu <- ppoints(1000)*5
postvar  <- 1/(n/sigsq + 1/tausqs)
postmean <- postvar*(mean(data)*n/sigsq + mu0/tausqs)
postdens <- matrix(NA,3,1000)
for (i in 1:3) {postdens[i,] <- dnorm(mu,postmean[i],sqrt(postvar[i]))}

plot(mu,postdens[1,],type="l",ylim=range(postdens),main="Posterior for Different Prior tausqs",col=1,lwd=2)
lines(mu,postdens[2,],col=2,lwd=2)
lines(mu,postdens[3,],col=3,lwd=2)
abline(v=mean(data),col=4,lwd=2)
abline(v=mu0,col=5,lwd=2)
legend(3.5,2.15,c("tausq=0.1","tausq=1","tausq=10","ybar","mu0"),col=1:5,lwd=2)

#ybar = 2.882
ybar = mean(data)
a=norm.post(1,0.10,ybar,100,5)
b=norm.post(1,1.00,ybar,100,5)
c=norm.post(1,10.0,ybar,100,5)


## repeat for a smaller sample size
n <- 10
data <- rnorm(n,mean=mu.true,sd=sqrt(sigsq))
par(mfrow=c(2,1),mar=c(3,3,2,1),las=1)
hist(data,main=paste("Normal Data (n=",length(data), ", ybar=",round(mean(data),3),")",sep=""))


mu <- ppoints(1000)*5
postvar  <- 1/(n/sigsq + 1/tausqs)
postmean <- postvar*(mean(data)*n/sigsq + mu0/tausqs)
postdens <- matrix(NA,3,1000)
for (i in 1:3) {postdens[i,] <- dnorm(mu,postmean[i],sqrt(postvar[i]))}

plot(mu,postdens[1,],type="l",ylim=range(postdens),col=1,main="Posterior for Different Prior tausqs",lwd=2)
lines(mu,postdens[2,],col=2,lwd=2)
lines(mu,postdens[3,],col=3,lwd=2)
abline(v=mean(data),col=4,lwd=2)
abline(v=mu0,col=5,lwd=2)
legend(3.5,1.35,c("tausq=0.1","tausq=1","tausq=10","ybar","mu0"),col=1:5,lwd=2)


#ybar = 2.399
ybar = mean(data)
a=norm.post(1,0.10,ybar,10,5)
b=norm.post(1,1.00,ybar,10,5)
c=norm.post(1,10.0,ybar,10,5)


## repeat for a larger sample size
n <- 1000
data <- rnorm(n,mean=mu.true,sd=sqrt(sigsq))
par(mfrow=c(2,1),mar=c(3,3,2,1),las=1)
hist(data,main=paste("Normal Data (n=",length(data), ", ybar=",round(mean(data),3),")",sep=""))


mu <- ppoints(1000)*5
postvar  <- 1/(n/sigsq + 1/tausqs)
postmean <- postvar*(mean(data)*n/sigsq + mu0/tausqs)
postdens <- matrix(NA,3,1000)
for (i in 1:3) {postdens[i,] <- dnorm(mu,postmean[i],sqrt(postvar[i]))}

plot(mu,postdens[1,],type="l",ylim=range(postdens),main="Posterior for Different Prior tausqs",col=1,lwd=2)
lines(mu,postdens[2,],col=2,lwd=2)
lines(mu,postdens[3,],col=3,lwd=2)
abline(v=mean(data),col=4,lwd=2)
abline(v=mu0,col=5,lwd=2)
legend(3.5,5.5,c("tausq=0.1","tausq=1","tausq=10","ybar","mu0"),col=1:5,lwd=2)


#ybar = 2.948
ybar = mean(data)
a=norm.post(1,0.10,ybar,1000,5)
b=norm.post(1,1.00,ybar,1000,5)
c=norm.post(1,10.0,ybar,1000,5)


#####################################################

## Reading in Data:
data <- read.table("players.post1970.txt",header=T,sep=",")
dim(data)

## Getting rid of pitchers
data <- data[data$pos != "P",]
dim(data)

## Reducing data to player-seasons where ab > 100
data <- data[data$AB > 100,]
dim(data)

## Getting rid of NA rows
problems <- which(is.na(data[,1]))
data <- data[-problems,]
dim(data)

## Calculating batting average
ba <- data$H/data$AB
n <- length(ba)
hist(ba)
min(ba)
data[which(ba==min(ba)),]
max(ba)
data[which(ba==max(ba)),]


## Function to sample from posterior for conjugate normal model


sample.norm.conj <- function(y,mu0,kappa0,nu0,sigsq0){
  n <- length(y)
  y.mean <- mean(y)
  y.ss <- var(y)*(n-1)
  discrep <- (y.mean-mu0)^2
  x <- rchisq(1,nu0+n)
  sigsq.samp <- (nu0*sigsq0+y.ss+kappa0*n*discrep/(kappa0+n))/x
  postvar <- 1/(n/sigsq.samp + kappa0/sigsq.samp)
  postmean <- postvar*(n*y.mean/sigsq.samp + kappa0*mu0/sigsq.samp)
  mu.samp <- rnorm(1,mean=postmean,sd=sqrt(postvar))
  out <- c(mu.samp,sigsq.samp)
  out
}

sample.norm.conj.many <- function(y,mu0,kappa0,nu0,sigsq0,numsamp){
  n <- length(y)
  y.mean <- mean(y)
  y.ss <- var(y)*(n-1)
  discrep <- (y.mean-mu0)^2
  x <- rchisq(numsamp,nu0+n)
  sigsq.samp <- (nu0*sigsq0+y.ss+kappa0*n*discrep/(kappa0+n))/x
  postvar <- 1/(n/sigsq.samp + kappa0/sigsq.samp)
  postmean <- postvar*(n*y.mean/sigsq.samp + kappa0*mu0/sigsq.samp)
  mu.samp <- rnorm(numsamp,mean=postmean,sd=sqrt(postvar))
  out <- cbind(mu.samp,sigsq.samp)
  out
}


## checking the posterior for different conjugate priors and non-informative 
##  (different values of mu0,kappa0)

theta1 <- sample.norm.conj.many(ba,0.0, 100,10,10,1000)  # mu0 = 0, kappa0 = 100, etc.
theta2 <- sample.norm.conj.many(ba,0.2, 100,10,10,1000)  # mu0 = 0.2, kappa0 = 100, etc.
theta3 <- sample.norm.conj.many(ba,0.0,1000,10,10,1000)  # mu0 = 0, kappa0 = 1000, etc.
theta4 <- sample.norm.conj.many(ba,0.2,   0, 0,10,1000)  # non-informative kappa0 = 0, nu0 = 0 


par(mfrow=c(4,1),mar=c(3,3,1,1))
xlim = range(theta1[,1],theta2[,1],theta3[,1],theta4[,1],mean(ba))
hist(theta1[,1],main="Mu: mu0=0,kappa0=100",xlim=xlim);   box(); abline(v=mean(ba),col=2)
hist(theta2[,1],main="Mu: mu0=0.2,kappa0=100",xlim=xlim); box(); abline(v=mean(ba),col=2)
hist(theta3[,1],main="Mu: mu0=0,kappa0=1000",xlim=xlim);  box(); abline(v=mean(ba),col=2)
hist(theta4[,1],main="Mu: non-informative",xlim=xlim);    box(); abline(v=mean(ba),col=2)



par(mfrow=c(4,1),mar=c(3,3,1,1))
xlim = range(theta1[,2],theta2[,2],theta3[,2],theta4[,2],var(ba))
hist(theta1[,2],main="Sigsq: mu0=0,kappa0=100",xlim=xlim);   box(); abline(v=var(ba),col=2)
hist(theta2[,2],main="Sigsq: mu0=0.2,kappa0=100",xlim=xlim); box(); abline(v=var(ba),col=2)
hist(theta3[,2],main="Sigsq: mu0=0,kappa0=1000",xlim=xlim);  box(); abline(v=var(ba),col=2)
hist(theta4[,2],main="Sigsq: non-informative",xlim=xlim);    box(); abline(v=var(ba),col=2)


## could have also generated mu directly from t distribution:
mu.sampt <- rt(10000,n-1)
mu.sampt <- mu.sampt*sqrt(var(ba)/n)+mean(ba)

## compare to original sampling scheme
theta <- sample.norm.conj.many(ba,0.2,0,0,10,10000)
mu.orig <- theta[,1]

par(mfrow=c(2,1))
xlim = range(mu.orig,mu.sampt)
hist(mu.orig,main="Mu: non-informative",xlim=xlim); abline(v=mean(ba),col=2)
hist(mu.sampt,main="Mu: from t dist",xlim=xlim);    abline(v=mean(ba),col=2)

