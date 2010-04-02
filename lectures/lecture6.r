#####################################################
########## STAT 289: CODE FOR LECTURE 6  ############
#####################################################


############################################
### GRID SEARCH FOR SIMPLE MIXTURE MODEL ###
############################################

## input data:  (mixture of normals data with true mus 0 and 2 and alpha=0.25)
data <- read.table("normnorm1.txt")
y <- data[,2]
n <- length(y)

## function for calculating log-likelihood:

loglik.mix <- function(alpha,mu1,mu0,y){
  phi1 <- dnorm(y,mu1,1)
  phi0 <- dnorm(y,mu0,1)
  loglik <- sum(log(alpha*phi1 + (1-alpha)*phi0))
  return(loglik)
}

loglike.mix <- function(theta,y){
  return(-loglik.mix(theta[1],theta[2],theta[3],y=y))
}


## calculating likelihood over a range of alpha, mu1, mu0
ngrid <- 21
alpha <- seq( 0,1,len=ngrid)   # alpha between 0 and 1
mu1 <- seq(-2,2,len=ngrid)   # mu1 between -2 and 2   
mu0 <- seq( 0,4,len=ngrid)   # mu0 between  0 and 4
loglike <- array(NA,dim=c(ngrid,ngrid,ngrid))

for (i in 1:ngrid){
  for (j in 1:ngrid){
    for (k in 1:ngrid){
      # impose the contraint mu1 > mu0 to make the model identifiable
      if(mu1[j]< mu0[k]) {loglike[i,j,k] <- loglik.mix(alpha[i],mu1[j],mu0[k],y)}
      if(mu1[j]>=mu0[k]) {loglike[i,j,k] <- -Inf}
    }
  }
}
like <- exp(loglike-max(loglike))


## plotting slices of likelihood as 2-d contour:
par(mfrow=c(2,2),las=1)
contour(mu1,mu0,loglike[ 4,,],main=paste("alpha=",alpha[ 4]),drawlabels=F)
contour(mu1,mu0,loglike[ 8,,],main=paste("alpha=",alpha[ 8]),drawlabels=F)
contour(mu1,mu0,loglike[12,,],main=paste("alpha=",alpha[12]),drawlabels=F)
contour(mu1,mu0,loglike[16,,],main=paste("alpha=",alpha[16]),drawlabels=F)


# Find the grid-based MLE 
theta <- as.matrix(expand.grid(alpha,mu1,mu0))
theta.grid <- theta[c(like)==max(c(like))]

newton <- nlminb(theta.grid,loglike.mix,y=y)
theta.mle <- newton$par

> theta.mle
[1] 0.3139937 0.1314180 2.0938594

> theta.grid
[1] 0.35 0.20 2.20






############################################
##### EM ALGORITHM FOR A MIXTURE MODEL #####
############################################

## input data: mixture of normals with true alpha=0.75,mu0=0,mu1=1,sigsq0=1,sigsq1=4
data <- read.table("normnorm2.txt")
y <- data[,2]
par(mfrow=c(1,1))
hist(y) 
n <- length(y)


## parameter names and true values
theta.label <- expression(alpha,mu[0],mu[1],sigma[0]^2,sigma[1]^2,loglike)
theta.true <- c(.75,0,1,1,4)

## function for calculating log-likelihood:

loglik.mixture <- function(alpha,mu0,mu1,sigsq0,sigsq1,y){
  phi0 <- dnorm(y,mu0,sqrt(sigsq0))
  phi1 <- dnorm(y,mu1,sqrt(sigsq1))
  loglik <- sum(log(alpha*phi1 + (1-alpha)*phi0))
  return(loglik)
}

loglike.mixture <- function(theta,y){ 
  return(-loglik.mixture(theta[1],theta[2],theta[3],theta[4],theta[5],y=y))
}


## Find MLE using Newton-Raphson
lower <- c(0,-1e10,-1e10,   0,   0)
upper <- c(1, 1e10, 1e10,1e10,1e10)
theta.start <- c(.5,-1,1,1,1)
newton <- nlminb(theta.true,loglike.mixture,y=y,lower=lower,upper=upper)
theta.nr <- c(newton$par,-newton$objective)




## Expectation function
Estep <- function(theta,y){
  n <- length(y)
  gamma <- rep(NA,n)
  for (i in 1:n){
    prob0 <- (1-theta[1])*dnorm(y[i],mean=theta[2],sd=sqrt(theta[4]))
    prob1 <-   theta[1]  *dnorm(y[i],mean=theta[3],sd=sqrt(theta[5]))
    gamma[i] <- prob1/(prob0+prob1)
  }
  return(gamma)
}

## Maximization function
Mstep <- function(gamma,y){
  n <- length(y)
  theta <- rep(NA,5)
  theta[1] <- sum(gamma)/n                                  # alpha
  theta[2] <- sum((1-gamma)*y)/sum(1-gamma)                 # mu0
  theta[3] <- sum(gamma*y)/sum(gamma)                       # mu1
  theta[4] <- sum((1-gamma)*((y-theta[2])^2))/sum(1-gamma)  # sigsq0
  theta[5] <- sum(gamma*((y-theta[3])^2))/sum(gamma)        # sigsq1
  return(theta)
}


## Starting values for EM algorithm:
theta <- c(0.5,-1,1,1,1)
loglike <- -loglike.mixture(theta,y)
itermat <- c(theta,loglike)

## Running EM iterations:
for (i in 1:100){
  gamma <- Estep(theta,y)
  theta <- Mstep(gamma,y)
  loglike <- -loglike.mixture(theta,y)
  itermat <- rbind(itermat,c(theta,loglike))
  cat(i,loglike,"\n")
}

## Looking at iterations
good <- 1:100
par(mfrow=c(3,2),mar=c(4,4,2,1),las=1)
for (i in c(6,1:5)){
  ylim <- range(itermat[good,i],theta.nr[i])
  plot(good,itermat[good,i],type="o",xlab="Iteration",ylab="",main=theta.label[i],ylim=ylim)
  abline(h=theta.nr[i],col=2,lwd=2)
}



## Better way to run EM
theta <- c(.5,-1,1,1,1)
loglike <- -loglike.mixture(theta,y)
itermat <- c(theta,loglike)
niter <- 1

repeat{
  niter <- niter + 1
  gamma <- Estep(theta,y)
  theta <- Mstep(gamma,y)
  loglike <- -loglike.mixture(theta,y)
  itermat <- rbind(itermat,c(theta,loglike))
  diff <- sum(abs(itermat[niter,1:5]-itermat[niter-1,1:5])/abs(itermat[niter,1:5])) 
  cat(niter,loglike,"\n")
  if (diff<1e-5 | niter>1000) break
}

> theta
> theta.nr[-6]
[1]  0.71700323 -0.06441126  0.98231351  1.20685213  4.23386001
[1]  0.71705675 -0.06443426  0.98224446  1.20664344  4.23377920

> -loglike.mixture(theta,y)
> -loglike.mixture(theta.nr,y)
[1] -2049.599
[1] -2049.599


## Looking at iterations
good <- 1:niter
par(mfrow=c(3,2),mar=c(4,4,2,1),las=1)
for (i in c(6,1:5)){
  ylim <- range(itermat[good,i],theta.nr[i])
  plot(good,itermat[good,i],type="o",xlab="Iteration",ylab="",main=theta.label[i],ylim=ylim)
  abline(h=theta.nr[i],col=2,lwd=2)
}





## Different Starting values: 
theta <- c(.5,-10,10,2,2)
loglike <- -loglike.mixture(theta,y)
itermat <- c(theta,loglike)
niter <- 1
repeat{
  niter <- niter + 1
  gamma <- Estep(theta,y)
  theta <- Mstep(gamma,y)
  loglike <- -loglike.mixture(theta,y)
  itermat <- rbind(itermat,c(theta,loglike))
  diff <- sum(abs(itermat[niter,1:5]-itermat[niter-1,1:5])/abs(itermat[niter,1:5])) 
  cat(niter,loglike,"\n")
  if (diff<1e-5 | niter>1000) break
}

## Looking at iterations
good <- 20:niter
par(mfrow=c(3,2),mar=c(4,4,2,1),las=1)
for (i in c(6,1:5)){
  ylim <- range(itermat[good,i],theta.nr[i])
  plot(good,itermat[good,i],type="o",xlab="Iteration",ylab="",main=theta.label[i],ylim=ylim)
  abline(h=theta.nr[i],col=2,lwd=2)
}



############################################
#### EM ALGORITHM FOR BASEBALL EXAMPLE #####
############################################

#Reading in Data:
data <- read.table("players.post1970.noP.txt",skip=1,sep=",",row.names=NULL)
hr <- data[,12]
ab <- data[,7]
year <- data[,2]
player <- data[,1]

#Reducing data to player-seasons where ab > 100
hr <- hr[ab>100]
year <- year[ab>100]
player <- player[ab>100]
ab <- ab[ab>100]

#Calculating homerun proportion:
hrprop <- hr/ab
par(mfrow=c(1,1))
hist(hrprop)


## Find the MLE using Newton-Raphson
lower <- c(0,-1e10,-1e10,   0,   0)
upper <- c(1, 1e10, 1e10,1e10,1e10)
theta.start <- c(.5,.03,.10,.01,.10)
newton <- nlminb(theta.start,loglike.mixture,y=hrprop,lower=lower,upper=upper)
theta.nr <- c(newton$par,-newton$objective)




#Running EM iterations
theta <- c(.1,.01,.07,.03,.03)
loglike <- -loglike.mixture(theta,hrprop)
itermat <- c(theta,loglike)

for (i in 1:200){
  gamma <- Estep(theta,hrprop)
  theta <- Mstep(gamma,hrprop)
  loglike <- -loglike.mixture(theta,hrprop)
  itermat <- rbind(itermat,c(theta,loglike))
  cat(i,loglike,"\n")
}



## Tracking iterations
theta.label <- expression(alpha,mu[0],mu[1],sigma[0]^2,sigma[1]^2,loglike)
good <- 1:201
par(mfrow=c(3,2),mar=c(4,4,2,1),las=1)
for (i in c(6,1:5)){
  ylim <- range(itermat[good,i],theta.nr[i])
  plot(good,itermat[good,i],type="o",xlab="Iteration",ylab="",main=theta.label[i],ylim=ylim)
  abline(h=theta.nr[i],col=2,lwd=2)
}


# plotting fitted mixture density
theta <- itermat[201,1:5]
par(mfrow=c(1,1))
hist(hrprop,prob=T)
x <- seq(0,.15,len=1001)
f0 <- (1-theta[1])*dnorm(x,theta[2],sqrt(theta[4]))
f1 <- theta[1] *dnorm(x,theta[3],sqrt(theta[5]))
lines(x,f0,col=2,lwd=2)
lines(x,f1,col=3,lwd=2)


# EM algorithm for equal-variance model
Estep2 <- function(theta,y){
  n <- length(y)  
  gamma <- rep(NA,n)
  for (i in 1:n){
    prob0 <- (1-theta[1])*dnorm(y[i],mean=theta[2],sd=sqrt(theta[4]))
    prob1 <-     theta[1]*dnorm(y[i],mean=theta[3],sd=sqrt(theta[4]))
    gamma[i] <- prob1/(prob0+prob1)
  }
  return(gamma)
}

Mstep2 <- function(gamma,y){
  n <- length(y)
  theta <- rep(NA,4)
  theta[1] <- sum(gamma)/n
  theta[2] <- sum((1-gamma)*y)/sum(1-gamma)
  theta[3] <- sum(gamma*y)/sum(gamma)
  theta[4] <- sum(gamma*(y-theta[3])^2 + (1-gamma)*(y-theta[2])^2)/n
  return(theta)
}

##observed data negative loglikelihood function for equal variance model
loglike.mixture2 <- function(theta,y){
  phi0 <- dnorm(y,theta[2],sqrt(theta[4]))
  phi1 <- dnorm(y,theta[3],sqrt(theta[4]))
  loglike <- sum(log(theta[1]*phi1 + (1-theta[1])*phi0))
  return(-loglike)
}


## Find the MLE using Newton-Raphson
lower <- c(0,-1e10,-1e10,   0)
upper <- c(1, 1e10, 1e10,1e10)
theta.start <- c(.2,.015,.06,.0002)
newton <- nlminb(theta.start,loglike.mixture2,y=hrprop,lower=lower,upper=upper)
theta.nr2 <- c(newton$par,-newton$objective)




#Running EM iterations
theta <- c(0.1,0.001,0.15,0.10)
loglike <- -loglike.mixture2(theta,hrprop)
itermat2 <- c(theta,loglike)

for (i in 1:200){
  gamma <- Estep2(theta,hrprop)
  theta <- Mstep2(gamma,hrprop)
  loglike <- -loglike.mixture2(theta,hrprop)
  itermat2 <- rbind(itermat2,c(theta,loglike))
  cat(i,loglike,"\n")
}


## Tracking iterations
theta.label <- expression(alpha,mu[0],mu[1],sigma^2,loglike)
good <- 10:201
par(mfrow=c(3,2),mar=c(4,4,2,1),las=1)
for (i in c(5,1:4)){
  ylim <- range(itermat2[good,i],theta.nr2[i])
  plot(good,itermat2[good,i],type="o",xlab="Iteration",ylab="",main=theta.label[i],ylim=ylim)
  abline(h=theta.nr2[i],col=2,lwd=2)
}


# plotting equal-variances fitted mixture density
theta <- itermat2[201,]
par(mfrow=c(1,1))
hist(hrprop,prob=T)
x <- seq(0,1,len=1001)
f0 <- (1-theta[1])*dnorm(x,theta[2],sqrt(theta[4]))
f1 <- theta[1]*dnorm(x,theta[3],sqrt(theta[4]))
lines(x,f0,col=2,lwd=2)
lines(x,f1,col=3,lwd=2)


#Getting Individual probabilities for each player

theta <- itermat2[201,]
gamma <- Estep2(theta,hrprop)

hist(gamma)
sum(gamma > 0.9999)
greatobs<-cbind(as.character(player[gamma>0.9999]),year[gamma>0.9999],hrprop[gamma>0.9999])



#####################################################
######## STAT 289: CODE FOR LECTURE 6, PART B #######
#####################################################


#####################################################
####### Data for Normal Hierarchical Model ##########
#####################################################

## simulated data with truesigsq <- 8, true tausq <- 1.5, true mu0 <- 0
data <- read.table("normhier.txt")
y<-data
par(mfrow=c(1,1))
boxplot(y) 

truesigsq <- 8 

##Calculating necessary statistics:
m <- length(y[1,])
n <- rep(NA,m)
means <- rep(NA,m)
for (i in 1:m){
  n[i] <- length(y[,i])
  means[i] <- mean(y[,i])
}
ntot <- sum(n)


#####################################################
# Sampling Parameters for Normal Hierarchical Model #
#####################################################

## finding right grid for tausq
tausq.grid <- ppoints(1000)*20

## function to evaluate marginal posterior for tausq
tausq.logpostfunc <- function(tausq){
  Vmu0 <- 1/sum(1/(tausq + truesigsq/n))
  mu0hat <- sum(means/(tausq + truesigsq/n))*Vmu0

  out <- 0.5*log(Vmu0)
  for (group in 1:m){
    out <- out - 0.5*log(tausq + truesigsq/n[group])
  }
  for (group in 1:m){
    out <- out - 0.5*((means[group]-mu0hat)^2)/(tausq + truesigsq/n[group])
  }
  out
}
tausq.logpost <- tausq.logpostfunc(tausq.grid)
tausq.post <- exp(tausq.logpost-max(tausq.logpost))
tausq.post <- tausq.post/sum(tausq.post)

par(mfrow=c(1,1))
plot(tausq.grid,tausq.post,type="l")

numsamp <- 1000
tausq.samp <- rep(NA,numsamp)
mu0.samp <- rep(NA,numsamp)
mu.samp <- matrix(NA,nrow=numsamp,ncol=m)
for (i in 1:numsamp){

  # sampling tausq from grid of values
  curtausq <- sample(tausq.grid,size=1,prob=tausq.post)

  # sampling mu0 given curtausq	
  Vmu0 <- 1/sum(1/(curtausq + truesigsq/n))
  mu0hat <- sum(means/(curtausq + truesigsq/n))*Vmu0
  curmu0 <- rnorm(1,mean=mu0hat,sd=sqrt(Vmu0))

  # sampling group means given curtausq and curmu0
  curmu <- rep(NA,m)
  for (j in 1:m){
    curvar <- 1/(n[j]/truesigsq + 1/curtausq)
    curmean <- (means[j]*n[j]/truesigsq + curmu0/curtausq)*curvar
    curmu[j] <- rnorm(1,mean=curmean,sd=sqrt(curvar))
  }
  tausq.samp[i] <- curtausq
  mu0.samp[i] <- curmu0
  mu.samp[i,] <- curmu
  print (i)
}

#####################################################
########### Examining Model Parameters ##############
#####################################################

par(mfrow=c(2,2))
hist(tausq.samp,main="tausq")
hist(mu0.samp,main="mu0")
hist(mu.samp[,1],main="mu group 1")
hist(mu.samp[,2],main="mu group 2")

# posterior probability group 5 has greater mean than group 1
postprob <- sum(mu.samp[,5] > mu.samp[,1])/numsamp

# posterior probability group 2 has greater mean than group 1
postprob <- sum(mu.samp[,2] > mu.samp[,1])/numsamp

#####################################################
######### Examining Shrinkage Graphically ###########
#####################################################

par(mfrow=c(1,1))
hist(y[,8],main="Eighth Group")
# add red line for posterior mean
abline(v=mean(mu.samp[,8]),col=2,lwd=2)
# add green line for group mean
abline(v=means[8],col=3,lwd=2)
# add blue line for overall mean
abline(v=mean(means),col=4,lwd=2)

#####################################################
######### Posterior Predictive Sampling #############
#####################################################

## sampling distribution of new observation 
## from a currently existing group

ystar.group1 <- rep(NA,numsamp)
ystar.group2 <- rep(NA,numsamp)
for (i in 1:numsamp){
  ystar.group1[i] <- rnorm(1,mean=mu.samp[i,1],sd=sqrt(truesigsq))
  ystar.group2[i] <- rnorm(1,mean=mu.samp[i,2],sd=sqrt(truesigsq))
}

par(mfrow=c(2,1))
xlim <- range(ystar.group1,ystar.group2)
hist(ystar.group1,main="Group 1 New Obs",xlim=xlim)
hist(ystar.group2,main="Group 2 New Obs",xlim=xlim)

## sampling distribution of new observation
## from an entirely new group

ystar.newgroup <- rep(NA,numsamp)
for (i in 1:numsamp){
  mu.newgroup <- rnorm(1,mean=mu0.samp[i],sd=sqrt(tausq.samp[i]))
  ystar.newgroup[i] <- rnorm(1,mean=mu.newgroup,sd=sqrt(truesigsq))
}

par(mfrow=c(3,1))
xlim <- range(ystar.group1,ystar.group2,ystar.newgroup)
hist(ystar.group1,main="Group 1 New Obs",xlim=xlim)
hist(ystar.group2,main="Group 2 New Obs",xlim=xlim)
hist(ystar.newgroup,main="New Group New Obs",xlim=xlim)

