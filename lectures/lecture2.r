############################################################
############## R CODE FOR STAT 289 LECTURE 2 ###############
############################################################

### Comparing Prior and Posterior for Theta_B ###
x <- ppoints(1000)
theta.b <- dbeta(x,27,23)
par(mfrow=c(1,1))
plot(x,theta.b,type="l",col=2,lwd=2,main="Distribution of Theta_B")
lines(x,rep(1,1000),type="l",col=1,lwd=2)
legend(0,5,c("prior","posterior"),col=c(1,2),lwd=2)


### Comparing Posterior for Theta_B and Posterior for Theta_W ###
x <- ppoints(1000)
theta.b <- dbeta(x,27,23)
theta.w <- dbeta(x,207,54)
plot(x,theta.w,type="l",col=4,lwd=2,main="Comparing Theta_B to Theta_W")
lines(x,theta.b,type="l",col=2,lwd=2)
lines(x,rep(1,1000),type="l",col=1,lwd=2)
legend(0,15,c("prior","theta.b","theta.w"),col=c(1,2,4),lwd=2)


### Simulating values to calculate posterior probability ###
theta.b.samp <- rbeta(10000,27,23)
theta.w.samp <- rbeta(10000,207,54)
plot(theta.b.samp,theta.w.samp,xlim=c(0,1),ylim=c(0,1))
abline(0,1)
postprob <- sum(theta.b.samp < theta.w.samp)/10000


### Simulating value from posterior predictive distribution ###
theta.b.samp <- rbeta(10000,27,23)
ystar.b.samp <- rbinom(10000,100,theta.b.samp)
hist(ystar.b.samp,main="Posterior Predictive Distribution of Ystar")


### Comparing posterior predictive distribution to Binomial based on MLE ###
theta.b.mle <- 26/48
ystar.alt.samp <- rbinom(10000,100,theta.b.mle)
par(mfrow=c(2,1))
hist(ystar.b.samp,main="Posterior Predictive Distribution of Ystar",xlim=c(0,100))
hist(ystar.alt.samp,main="Parametric Bootstrap Distribution of Ystar",xlim=c(0,100))


##################################################################

### Comparing Different Priors for Binomial probability theta ###
theta <- ppoints(1000)
density.lapl <- dbeta(theta,1,1)
density.jeff <- dbeta(theta,0.5,0.5)
density.mle <- dbeta(theta,0,0)
density.mle <- dbeta(theta,0.000001,0.000001)
ylim = range(density.lapl,density.jeff,density.mle)
par(mfrow=c(1,1))
plot(theta,density.lapl,type="l",col=1,lwd=2,main="Prior Distributions of Theta",ylim=ylim)
lines(theta,density.jeff,type="l",col=2,lwd=2)
lines(theta,density.mle,type="l",col=3,lwd=2)
legend(0.4,10,c("Beta(1,1)","Beta(0.5,0.5)","Beta(0,0)"),col=c(1,2,3),lwd=2)


### Comparing Posterior for Different Priors for Small Dataset ###
### Small Dataset: 4 successes, 2 failures
theta <- ppoints(1000)
posterior.lapl <- dbeta(theta,5,3)
posterior.jeff <- dbeta(theta,4.5,2.5)
posterior.mle <- dbeta(theta,4.000001,2.000001)
ylim = range(posterior.lapl,posterior.jeff,posterior.mle)
plot(theta,posterior.lapl,type="l",col=1,lwd=2,main="Posterior Distributions: Small Dataset",ylim=ylim)
lines(theta,posterior.jeff,type="l",col=2,lwd=2)
lines(theta,posterior.mle,type="l",col=3,lwd=2)
legend(0,2.3,c("Prior=Beta(1,1)","Prior=Beta(0.5,0.5)","Prior=Beta(0,0)"),col=c(1,2,3),lwd=2)
 
### Comparing Posterior distributions to Sampling distribution using CLT ###
theta.mle <- 4/6
theta.mle.var <- theta.mle*(1-theta.mle)/6
theta.mle.sd <- sqrt(theta.mle.var)
theta.sampdist <- dnorm(theta,mean=theta.mle,sd=theta.mle.sd)
plot(theta,posterior.lapl,type="l",col=1,lwd=2,main="Posterior Distributions: Small Dataset",ylim=ylim)
lines(theta,posterior.jeff,type="l",col=2,lwd=2)
lines(theta,posterior.mle,type="l",col=3,lwd=2)
lines(theta,theta.sampdist,type="l",col=4,lwd=2)
legend(0,2.3,c("Prior=Beta(1,1)","Prior=Beta(0.5,0.5)","Prior=Beta(0,0)","Sampling Dist"),col=c(1,2,3,4),lwd=2)


### Comparing Posterior for Different Priors for Teal Blacks Data ###
### Remember Blacks in Teal example: 26 successes, 22 failures
theta <- ppoints(1000)
posterior.lapl <- dbeta(theta,27,23)
posterior.jeff <- dbeta(theta,26.5,22.5)
posterior.mle <- dbeta(theta,26.000001,22.000001)
ylim = range(posterior.lapl,posterior.jeff,posterior.mle)
plot(theta,posterior.lapl,type="l",col=1,lwd=2,main="Posterior Distributions: Teal Blacks",ylim=ylim)
lines(theta,posterior.jeff,type="l",col=2,lwd=2)
lines(theta,posterior.mle,type="l",col=3,lwd=2)
legend(-.025,5.75,c("Prior=Beta(1,1)","Prior=Beta(0.5,0.5)","Prior=Beta(0,0)"),col=c(1,2,3),lwd=2)


### Comparing Posterior distributions to Sampling distribution using CLT ###
theta.mle <- 26/48
theta.mle.var <- theta.mle*(1-theta.mle)/48
theta.mle.sd <- sqrt(theta.mle.var)
theta.sampdist <- dnorm(theta,mean=theta.mle,sd=theta.mle.sd)
plot(theta,posterior.lapl,type="l",col=1,lwd=2,main="Posterior Distributions: Teal Blacks",ylim=ylim)
lines(theta,posterior.jeff,type="l",col=2,lwd=2)
lines(theta,posterior.mle,type="l",col=3,lwd=2)
lines(theta,theta.sampdist,type="l",col=4,lwd=2)
legend(-.025,5.75,c("Prior=Beta(1,1)","Prior=Beta(0.5,0.5)","Prior=Beta(0,0)","Sampling Dist"),col=c(1,2,3,4),lwd=2)




### Analysis for Nonconjugate Prior with Small Dataset (4 successes, 2 failures) ###
theta = ppoints(1000)
prior.tri = 2-4*abs(.5-theta)
likelihood = dbinom(4,6,theta)
posterior.tri = prior.tri*likelihood
theta.draws = sample(theta,10000,replace=T,prob=posterior.tri)
pdf("lecture04-plot6.pdf",height=8,width=6)
par(mfrow=c(4,1),mar=c(3,3,2,1))
plot(theta,prior.tri,type="l",main="Prior",lwd=2,col=4)
plot(theta,likelihood,type="l",main="Likelihood",lwd=2,col=4)
plot(theta,posterior.tri,type="l",main="Posterior",lwd=2,col=4)
hist(theta.draws,nclass=100,main="Posterior Samples",xlim=c(0,1))
dev.off()




### Ploting Posterior for Nonconjugate Prior for Teal Blacks Data ###
### Remember Blacks in Teal example: 26 successes, 22 failures
theta = ppoints(1000)
prior.tri = 2-4*abs(.5-theta)
likelihood = dbinom(26,48,theta)
posterior.tri = prior.tri*likelihood
posterior.tri = posterior.tri/sum(posterior.tri)
pdf("lecture04-plot7.pdf",height=7,width=6)
par(mfrow=c(3,1),mar=c(3,3,2,1))
plot(theta,prior.tri,type="l",main="Prior")
plot(theta,likelihood,type="l",main="Likelihood")
plot(theta,posterior.tri,type="l",main="Posterior")
dev.off()


moments.beta = function(alpha,beta){
  mean = alpha/(alpha+beta)
  sdev = sqrt(alpha*beta/((alpha+beta)^2*(alpha+beta+1)))
  mode = (alpha-1)/(alpha+beta-2)
  return(list(mean=mean,sdev=sdev,mode=mode))
}

moments.beta(1,1)
moments.beta(.5,.5)
moments.beta(0,0)

moments.beta(1+4,1+2)
moments.beta(.5+4,.5+2)
moments.beta(0+4,0+2)

moments.beta(1+26,1+22)
moments.beta(.5+26,.5+22)
moments.beta(0+26,0+22)



##### simulating some Normal data ####
mu.true <- 3
sigsq <- 5
n <- 100
data <- rnorm(n,mean=mu.true,sd=sqrt(sigsq))
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



## repeat for a smaller sample size
n <- 10
data <- rnorm(n,mean=mu.true,sd=sqrt(sigsq))
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



## repeat for a larger sample size
n <- 1000
data <- rnorm(n,mean=mu.true,sd=sqrt(sigsq))
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


norm.post(1,0.1,2.882,100,5)
norm.post(1,1.0,2.882,100,5)
norm.post(1,10.,2.882,100,5)

norm.post(1,0.1,2.399,10,5)
norm.post(1,1.0,2.399,10,5)
norm.post(1,10.,2.399,10,5)

norm.post(1,0.1,2.948,1000,5)
norm.post(1,1.0,2.948,1000,5)
norm.post(1,10.,2.948,1000,5)
