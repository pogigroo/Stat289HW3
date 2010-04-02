#####################################################
########## STAT 289: CODE FOR LECTURE 7  ############
#####################################################

#####################################################
## Binomial Hierarchical Model with Rat Tumor Data ##
#####################################################

data <- read.table("rattumors.txt",header=T)
y <- data[,1]
n <- data[,2]
m <- length(y)
hist(y/n)


logpost.ab <- function(alpha,beta,y,n){
  ######################################################
  # Return the log marginal posterior for (alpha,beta)
  ######################################################
  logprior <- -5/2*log(alpha+beta)
  loglike <- sum(lbeta(alpha+y,beta+n-y)-lbeta(alpha,beta))
  logpost <- logprior + loglike
  return(logpost)
}


alphagrid <- ppoints(100)*5  # alpha between 0 and 5
betagrid <- ppoints(100)*30  # beta between 0 and 30
logpost <- matrix(NA,100,100)
for (i in 1:100){ 
  for (j in 1:100){
    logpost[i,j] <- logpost.ab(alphagrid[i],betagrid[j],y,n)
  }
}
post <- exp(logpost-max(logpost))
post <- post/sum(post)
contour(alphagrid,betagrid,post,drawlabels=F)

post.alpha <- apply(post,1,sum)
plot(alphagrid,post.alpha,xlab=expression(alpha),type="l",col=4)

post.betacond <- matrix(NA,100,100)
for (i in 1:100) post.betacond[i,] <- post[i,]/post.alpha[i]

alphagrid[25]
alphagrid[50]
alphagrid[75]
par(mfrow=c(3,1))
plot(betagrid,post.betacond[25,],type="l",main="dist. of beta for alpha = 1.225")
plot(betagrid,post.betacond[50,],type="l",main="dist. of beta for alpha = 2.475")
plot(betagrid,post.betacond[75,],type="l",main="dist. of beta for alpha = 3.725")

numsamp <- 1000
alpha.samp <- rep(NA,numsamp)
beta.samp <- rep(NA,numsamp)
theta.samp <- matrix(NA,nrow=numsamp,ncol=m)
for (iter in 1:numsamp){
  a <- sample(1:length(alphagrid),size=1,prob=post.alpha)
  b <- sample(1:length(betagrid),size=1,prob=post.betacond[a,])
  alpha.samp[iter] <- alphagrid[a]
  beta.samp[iter] <- betagrid[b]
  theta.samp[iter,] <- rbeta(m,alpha.samp[iter]+y,beta.samp[iter]+n-y)
}

par(mfrow=c(3,1))
hist(alpha.samp)
hist(beta.samp)
hist(alpha.samp/(alpha.samp+beta.samp))

### Examining Shrinkage

theta.postmean <- apply(theta.samp,2,mean)
overall.postmean <- mean(alpha.samp/(alpha.samp+beta.samp))

par(mfrow=c(1,1))
plot(1:71,y/n,main="Shrinkage of Sample Proportions",pch=19)
abline(h=overall.postmean,col=4,lwd=2)
points(1:71,theta.postmean,pch=19,col=2)
legend(2,0.35,c("Sample Prop","Post Mean","Overall Mean"),pch=19,col=c(1,2,4))


#####################################################
####### Data for Normal Hierarchical Model ##########
#####################################################

## simulated data with truesigsq = 8, true tausq = 1.5, true mu0 = 0
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
######## EM for Normal Hierarchical Model ###########
#####################################################

itermat <- matrix(NA,nrow=1,ncol=2)
itermat[1,1] <- -10  # starting value for mu0 
itermat[1,2] <- 5    # starting value for tausq
diff <- 1
numiters <- 1
while (diff > 0.000001){
  mu0 <- itermat[numiters,1]
  tausq <- itermat[numiters,2]

  # Expection step
  temp <- n/truesigsq + 1/tausq 
  Emu <- (n*means/truesigsq + mu0/tausq)/temp
  Emu2 <- Emu^2 + 1/temp

  # Maximization function:
  mu0new <- sum(Emu)/m
  tausq0new <- sum(Emu2-2*mu0*Emu+mu0^2)/m
  itermat <- rbind(itermat,c(mu0new,tausq0new))
  numiters <- numiters + 1
  print (numiters)
  diff <- sum(abs(itermat[numiters,]-itermat[numiters-1,])/abs(itermat[numiters,])) 
}

mu0.max <- itermat[numiters,1]
tausq.max <- itermat[numiters,2]

temp <- n/truesigsq + 1/tausq.max
muvec.max <- (n*means/truesigsq + mu0.max/tausq.max)/temp

par(mfrow=c(1,1))
plot(1:10,means,main="Shrinkage of Normal Means",pch=19)
abline(h=mu0.max,col=4,lwd=2)
points(1:10,muvec.max,pch=19,col=2)
legend(8,2,c("Data Mean","Post Mean","Mu0"),pch=19,col=c(1,2,4))



#####################################################
####### STAT 289: CODE FOR LECTURE 7, part b ########
#####################################################

#################################################
########## Simple Random Walk ###################
#################################################

theta <- rep(NA,20)
theta[1] <- 0
for (i in 2:20){
  u <- runif(1,0,1)
  if (u <= 0.5) {theta[i] <- theta[i-1]+1}
  if (u > 0.5)  {theta[i] <- theta[i-1]-1}
}
par(mfrow=c(1,1))
plot(1:20,theta,type="l",col=1,ylim=c(-10,10))


#################################################
############# Brownian Motion ###################
#################################################

theta <- matrix(NA,1000,8)
for (j in 1:8){
  theta[1,j] <- 0
  for (i in 2:1000) theta[i,j] <- rnorm(1,theta[i-1,j],1)
}
matplot(theta,type="l")


#################################################
############### Stationary AR(1) ################
#################################################

T <- 1000
theta <- matrix(NA,T,3)
theta.start <- c(-20,0,20)
for (j in 1:3){
  theta[1,j] <- theta.start[j]
  for (i in 2:T) theta[i,j] <- rnorm(1,.9*theta[i-1,j],1)
}
matplot(theta,type="l")
points(rep(1,3),theta.start,col=1:3,pch=16,cex=2)
 


#################################################
###### BIVARIATE NORMAL GIBBS SAMPLER ###########
#################################################


## parameter values:
mu1 <- 0
mu2 <- 0
sigsq1 <- 1
sigsq2 <- 1
rho <- 0.8

## initial values:
theta1 <- rep(NA,1000)
theta2 <- rep(NA,1000)
theta1[1] <- 10
theta2[1] <- 10

## Gibbs sampler:
for (i in 2:1000){
  theta1[i] <- rnorm(1,rho*theta2[i-1],sqrt(1-rho^2))
  theta2[i] <- rnorm(1,rho*theta1[i],sqrt(1-rho^2))
}



## checking univariate means and variances:
means1 <- rep(NA,1000)
means2 <- rep(NA,1000)
vars1 <- rep(NA,1000)
vars2 <- rep(NA,1000)
for (i in 1:1000){
  means1[i] <- mean(theta1[1:i])
  means2[i] <- mean(theta2[1:i])
  vars1[i] <- var(theta1[1:i])
  vars2[i] <- var(theta2[1:i])
}

par(mfrow=c(2,2))
plot(1:1000,means1,type="l",main="mean(theta1)")
abline(h=0,col=2)
plot(1:1000,means2,type="l",main="mean(theta2)")
abline(h=0,col=2)
plot(1:1000,vars1,type="l",main="var(theta1)")
abline(h=1,col=2)
plot(1:1000,vars2,type="l",main="var(theta2)")
abline(h=1,col=2)

## checking convergence by comparing multiple chains from different start points:
theta1.alt <- rep(NA,1000)
theta2.alt <- rep(NA,1000)
theta1.alt[1] <- -10
theta2.alt[1] <- -10
for (i in 2:1000){
  theta1.alt[i] <- rnorm(1,rho*theta2.alt[i-1],sqrt(1-rho^2))
  theta2.alt[i] <- rnorm(1,rho*theta1.alt[i],sqrt(1-rho^2))
}

par(mfrow=c(2,1))
plot(1:1000,theta1,col=1,ylim=range(theta1,theta1.alt),type="l")
lines(1:1000,theta1.alt,col=2)
plot(1:1000,theta2,col=1,ylim=range(theta2,theta2.alt),type="l")
lines(1:1000,theta2.alt,col=2)

## throwing out "burn-in" samples:

theta1.new <- theta1[-c(1:100)]
theta2.new <- theta2[-c(1:100)]

## plotting auto-correlation of each parameter

par(mfrow=c(2,1))
acf(theta1.new)
acf(theta2.new) 

## retaining independent samples by only taking every fifth sample:

temp <- 5*c(1:(900/5))
theta1.final <- theta1.new[temp]
theta2.final <- theta2.new[temp]

par(mfrow=c(2,1))
acf(theta1.final)
acf(theta2.final) 



plot.bivariate.normal.draws <- function(theta1,theta2,iter,rho,newplot=T,col=4){
  ################################################################
  # Plots results of Gibbs sampler for bivariate normal
  # 1. Contour plot of bivariate normal 
  # 2. Superimpose draws from Gibbs Sampler
  ################################################################
  if (newplot){
    dmnorm <- function(x,m,V){return(det(V)^(-.5)*exp(-.5*t(x-m)%*%solve(V)%*%(x-m)))}

    x1 <- x2 <- seq(-4,4,len=50)
    pdf <- matrix(NA,50,50)

    for (i in 1:50) 
      for (j in 1:50) 
        pdf[i,j] <- dmnorm(c(x1[i],x2[j]),m=c(0,0),V=matrix(c(1,rho,rho,1),2,2))

    par(mfrow=c(1,1),las=1)
    contour(x1,x2,pdf)
  }
  lines(theta1[iter],theta2[iter],col=col)
}

## Bivariate plot of Gibbs sampler draws
plot.bivariate.normal.draws(theta1,theta2,1:1,rho,newplot=T)
plot.bivariate.normal.draws(theta1,theta2,1:10,rho,newplot=F)
plot.bivariate.normal.draws(theta1,theta2,1:100,rho,newplot=F)
plot.bivariate.normal.draws(theta1,theta2,1:1000,rho,newplot=F)


#########################################
## What about a correlation of 0.99 ?  ##
#########################################

rho <- 0.99

theta1 <- rep(NA,1000)
theta2 <- rep(NA,1000)
theta1[1] <- 10
theta2[1] <- 10
theta1.alt <- rep(NA,1000)
theta2.alt <- rep(NA,1000)
theta1.alt[1] <- -10
theta2.alt[1] <- -10

## Gibbs sampler:
for (i in 2:1000){
  theta1[i] <- rnorm(1,rho*theta2[i-1],sqrt(1-rho^2))
  theta2[i] <- rnorm(1,rho*theta1[i],sqrt(1-rho^2))
}
for (i in 2:1000){
  theta1.alt[i] <- rnorm(1,rho*theta2.alt[i-1],sqrt(1-rho^2))
  theta2.alt[i] <- rnorm(1,rho*theta1.alt[i],sqrt(1-rho^2))
}

par(mfrow=c(2,1))
plot(1:1000,theta1,col=1,ylim=range(theta1,theta1.alt),type="l")
lines(1:1000,theta1.alt,col=2)
plot(1:1000,theta2,col=1,ylim=range(theta2,theta2.alt),type="l")
lines(1:1000,theta2.alt,col=2)


## Bivariate plot of Gibbs sampler draws
plot.bivariate.normal.draws(theta1,theta2,1:1,rho,newplot=T)
plot.bivariate.normal.draws(theta1,theta2,1:100,rho,newplot=F)
plot.bivariate.normal.draws(theta1,theta2,1:200,rho,newplot=F)
plot.bivariate.normal.draws(theta1,theta2,1:1000,rho,newplot=F)


## throwing out "burn-in" samples:

theta1.new <- theta1[-c(1:200)]
theta2.new <- theta2[-c(1:200)]

## plotting auto-correlation of each parameter

par(mfrow=c(2,1))
acf(theta1.new)
acf(theta2.new) 

