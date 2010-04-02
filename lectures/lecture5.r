##############################################################
############## STAT 289: CODE FOR LECTURE 5  #################
##############################################################

par(mfrow=c(1,1),mar=c(3.5,3.5,2,1),mgp=c(2,.65,0),las=1)

##############################################################
########## UNIVARIATE NEWTON-RAPHSON: GAMMA DATA #############
##############################################################

## input data:  (gamma data with true alpha <- 2 and beta <- 3)
y <- read.table("gammaexample.txt")
y <- y[,2]
n <- length(y)


## Negative loglikelihood, gradient, hessian for (alpha,beta)
loglike.gamma <- function(theta,y){
  n <- length(y)
  f <- n*theta[1]*log(theta[2]) - n*lgamma(theta[1])
  f <- f + (theta[1]-1)*sum(log(y)) - theta[2]*sum(y)
  return(-f)
}
gradient.gamma <- function(theta,y){
  n <- length(y)
  g <- rep(NA,2)
  g[1] <- n*log(theta[2]) - n*digamma(theta[1]) + sum(log(y))
  g[2] <- n*theta[1]/theta[2] - sum(y)
  return(-g)
}
hessian.gamma <- function(theta,y){
  n <- length(y)
  H <- matrix(NA,2,2)
  H[1,1] <- -n*trigamma(theta[1])
  H[1,2] <- n/theta[2]
  H[2,1] <- n/theta[2]
  H[2,2] <- -n*theta[1]/theta[2]^2
  return(-H)
}
log.like.gamma <- function(alpha,beta,y){
  n <- length(y)
  value <- n*alpha*log(beta) - n*lgamma(alpha)
  value <- value + (alpha-1)*sum(log(y)) - beta*sum(y)
  return(value)
}
loglike.alpha <- function(alpha,beta,y){
  return(-log.like.gamma(alpha,beta,y))
}
gradient.alpha <- function(alpha,beta,y){
  return(-(length(y)*(log(beta)-digamma(alpha))+sum(log(y))))
}
hessian.alpha <- function(alpha,beta,y){
  return(as.matrix(length(y)*trigamma(alpha)))
}


contour.likelihood.gamma <- function(amin,amax,bmin,bmax,ngrid,y){
  ############################################################
  ## Contour plot of likelihood for gamma sampling model
  ############################################################
  alpha <- seq(amin,amax,len=ngrid)
  beta  <- seq(bmin,bmax,len=ngrid)
  loglike <- matrix(NA,ngrid,ngrid)
  for (i in 1:ngrid){
    for (j in 1:ngrid){
      loglike[i,j] <- log.like.gamma(alpha[i],beta[j],y)
    }
  }
  like <- exp(loglike - max(loglike))
  main <- expression("Likelihood function for "*alpha*" and "*beta*"")
  contour(alpha,beta,like,drawlabels=F,xlab=expression(alpha),ylab=expression(beta),main=main)
  mle <- expand.grid(alpha,beta)[c(like)==max(c(like)),]
  return(list(alpha=alpha,beta=beta,like=like,mle=mle))
}

plot.likelihood.alpha <- function(amin,amax,beta,ngrid,y){
  #######################################################
  # Plot likelihood of alpha for fixed beta
  #######################################################
  alpha <- seq(amin,amax,len=ngrid)
  loglike <- rep(NA,ngrid)
  for (i in 1:ngrid){
    loglike[i] <- log.like.gamma(alpha[i],beta,y)
  }
  like <- exp(loglike - max(loglike))
  main <- substitute(paste("Likelihood for "*alpha*" when "*beta*"=",b),list(b=beta))
  plot(alpha,like,type="l",xlab=expression(alpha),ylab="",main=main)
  mle <- alpha[like==max(like)]
  return(list(alpha=alpha,like=like,mle=mle))
}

plot.gradient.alpha <- function(amin,amax,beta,ngrid,y){
  #######################################################
  # Plot gradient of minus loglikelihood of alpha 
  #######################################################
  alpha <- seq(amin,amax,len=ngrid)
  gradient <- rep(NA,ngrid)
  for (i in 1:ngrid){
    gradient[i] <- gradient.alpha(alpha[i],beta,y)
  }
  main <- substitute(paste("-Score Function for "*alpha*" when "*beta*"=",b),list(b=beta))
  plot(alpha,gradient,type="l",xlab=expression(alpha),ylab="",main=main)
  return(list(alpha=alpha,gradient=gradient))
}

newton.raphson.alpha <- function(alpha0,beta,y,abs.tol=1e-6,max.iter=10){
  #######################################################
  # Run a Newton-Raphson algorithm for alpha
  #######################################################
  iter <- 1
  alpha <- alpha0
  f <- loglike.alpha(alpha,beta=beta,y=y)
  alpha.iter <- alpha
  repeat{
    iter <- iter+1
    g <- gradient.alpha(alpha,beta=beta,y=y)
    H <- hessian.alpha(alpha,beta=beta,y=y)
    alpha.new <- alpha - solve(H)%*%g
    if (alpha.new<=0) alpha.new=1e-2
    f.new <- loglike.alpha(alpha.new,beta=beta,y=y)
    f.diff <- abs(f.new - f)
    f <- f.new
    alpha <- alpha.new
    alpha.iter <- c(alpha.iter,alpha)
    if (f.diff<abs.tol | iter>=max.iter) break
  }
  return(list(alpha=alpha.iter,mle=alpha))
}

## Plot posterior (=likelihood) for (alpha,beta)
par(mfrow=c(1,1))
grid <- contour.likelihood.gamma(1,3,1,5,100,y)
grid <- contour.likelihood.gamma(1.5,4.5,2,6,100,y)
grid <- contour.likelihood.gamma(1.5,3.5,2,5.5,100,y)
points(grid$mle[1],grid$mle[2],pch="X",col=4,cex=2)

## Plot posterior for alpha for beta=3
par(mfrow=c(2,1))
par(ask=F)
grid <- contour.likelihood.gamma(1.5,3.5,2,5.5,100,y)
abline(3,0,lty=2,col=4)
a <- plot.likelihood.alpha(1.5,3.5,100,y=y,beta=3)
a <- plot.likelihood.alpha(1.5,2.5,100,y=y,beta=3)
b <- plot.gradient.alpha(1.5,2.5,100,y=y,beta=3)
abline(h=0,col=4,lty=2)

## Use built-in R function nlminb to find the mode
alpha0 <- 10
newton1 <- nlminb(alpha0,loglike.alpha,beta=3,y=y,lower=1e-5)
newton2 <- nlminb(alpha0,loglike.alpha,gradient.alpha,hessian.alpha,beta=3,y=y,lower=1e-5)
abline(v=newton1$par,col=2)
abline(v=newton2$par,col=4)

## Posterior mode and asymptotic posterior intervals
alpha.hat <- newton1$par
var.hat <- solve(hessian.alpha(alpha.hat,beta=3,y))
alpha.est <- alpha.hat[1] + c(-1,0,1)*1.96*sqrt(diag(var.hat)[1])
alpha.est


## Try the Newton-Raphson from four different starting values 
par(mfrow=c(1,1))
par(ask=F)
alpha.start <- c(.1,3,5,10)
plot(1,1,type="n",xlab="Iteration",ylab=expression(alpha),xlim=c(1,12),
     ylim=range(alpha.start),main=expression("Starting Values for "*alpha*""))
points(rep(1,4),alpha.start,col=2,pch=16,cex=2)
abline(h=newton1$par,col=1)


## Track iterations of newton-raphson from different starting values
## (function newton.raphson.alpha is defined above)
par(ask=T)
for (i in 1:4){
  plot(1,1,type="n",xlab="Iteration",ylab=expression(alpha),xlim=c(1,12),
       ylim=range(alpha.start),main=expression("Newton-Raphson Iterations for "*alpha*""))
  points(1,alpha.start[i],col=2,pch=16,cex=2)
  abline(h=newton1$par,col=1)
  newton <- newton.raphson.alpha(alpha.start[i],beta=3,y,max.iter=12)
  lines(newton$alpha,type="o",col=2,pch=16)
}



##############################################################
######### MULTIVARIATE NEWTON-RAPHSON: GAMMA DATA ############
##############################################################

# Plot posterior and show the grid-based posterior mode
par(ask=F)
par(mfrow=c(1,1))
grid <- contour.likelihood.gamma(1,3,1,5,100,y)
grid <- contour.likelihood.gamma(1.5,4.5,2,6,100,y)
grid <- contour.likelihood.gamma(1.5,3.5,2,5.5,100,y)
points(grid$mle[1],grid$mle[2],pch="X",col=4,cex=2)

## Use built-in R function nlminb to find the mode
theta0 <- c(3,2.5)
newton1 <- nlminb(theta0,loglike.gamma,y=y,lower=rep(1e-5,2))
newton2 <- nlminb(theta0,loglike.gamma,gradient.gamma,hessian.gamma,y=y,lower=rep(1e-5,2))

# Plot the starting values and posterior mode
grid <- contour.likelihood.gamma(1.5,3.5,2,5.5,100,y)
points(theta0[1],theta0[2],pch=16,col=2,cex=2)
points(newton1$par[1],newton1$par[2],pch="X",col=2,cex=2)
points(newton2$par[1],newton2$par[2],pch="X",col=4,cex=2)


# Posterior mode and asymptotic posterior intervals
theta.hat <- newton1$par
covar.hat <- solve(hessian.gamma(theta.hat,y))
alpha.est <- theta.hat[1] + c(-1,0,1)*1.96*sqrt(diag(covar.hat)[1])
beta.est  <- theta.hat[2] + c(-1,0,1)*1.96*sqrt(diag(covar.hat)[2])
alpha.est
beta.est


newton.raphson.gamma <- function(theta0,y,abs.tol=1e-6,max.iter=10){
  #################################################################  
  # Run a Newton-Raphson algorithm in the gamma example
  #################################################################
  iter <- 1
  theta <- theta0
  f <- loglike.gamma(theta,y)
  theta.iter <- NULL
  theta.iter <- rbind(theta.iter,t(theta))
  repeat{
    iter <- iter+1
    g <- gradient.gamma(theta,y)
    H <- hessian.gamma(theta,y)
    theta.new <- theta - solve(H)%*%g
    if (theta.new[1]<=0) theta.new[1]=1e-2
    if (theta.new[2]<=0) theta.new[2]=1e-2
    f.new <- loglike.gamma(theta.new,y)
    f.diff <- abs(f.new - f)
    f <- f.new
    theta <- theta.new
    theta.iter <- rbind(theta.iter,t(theta))
    if (f.diff<abs.tol | iter>=max.iter) break
  }
  return(list(theta=theta.iter,mle=theta))
}


# Try newton-raphson from four different starting values 
par(ask=F)
theta <- rbind(c(1.5,2),c(3.5,2),c(3.5,5.5),c(1.5,5.5))
grid <- contour.likelihood.gamma(0.1,4,0.1,5.5,100,y)
points(theta[,1],theta[,2],col=2,pch=16,cex=2)


## Track newton-raphson interations for different starting values
## (function newton.raphson.gamma is defined above)
par(ask=T)
for (i in 1:4){
  newton <- newton.raphson.gamma(theta[i,],y,max.iter=15)
  grid <- contour.likelihood.gamma(0.1,4,0.1,5.5,100,y)
  points(theta[i,1],theta[i,2],col=2,pch=16,cex=2)
  lines(newton$theta[,1],newton$theta[,2],type="o",col=2,pch=16)
}



################################################################
## Example 2: Multivariate Newton-Raphson: Poisson Planes 
################################################################

## input data:  
data <- read.table("planes.txt",skip=1)
accid <- data[,2]
time <- data[,1]-1976

n <- length(accid)
#pdf("lecture08-plot6.pdf",height=5,width=6)
plot(time,accid)
lm <- lm(accid~time)
abline(lm)
#dev.off()


################################################################
## Negative loglikelihood, gradient, hessian for planes data
################################################################

log.like.planes <- function(alpha,beta,y,t){
  if (alpha+beta*max(t)<=0) value <- -Inf
  if (alpha+beta*max(t)>0) value <- sum(y*log(alpha+beta*t)-(alpha+beta*t))
  return(value)
}
loglike.planes <- function(theta,y,t){
  return(-log.like.planes(theta[1],theta[2],y=y,t=t))
}
gradient.planes <- function(theta,y,t){
  x <- theta[1]+theta[2]*t
  g <- rep(NA,2)
  g[1] <- sum(y/x) - length(y)
  g[2] <- sum(y*t/x) - sum(t)
  return(-g)
}
hessian.planes <- function(theta,y,t){
  x <- theta[1]+theta[2]*t
  H <- matrix(NA,2,2)
  H[1,1] <- -sum(y/x^2)
  H[1,2] <- -sum(y*t/x^2)
  H[2,1] <- -sum(y*t/x^2)
  H[2,2] <- -sum(y*t^2/x^2)
  return(-H)
}


contour.likelihood.planes <- function(amin,amax,bmin,bmax,ngrid,y,t){
  #################################################################  
  # Contour plot of likelihood in poisson planes model
  #################################################################
  alpha <- seq(amin,amax,len=ngrid)
  beta  <- seq(bmin,bmax,len=ngrid)
  loglike <- matrix(NA,ngrid,ngrid)
  for (i in 1:ngrid){
    for (j in 1:ngrid){
      loglike[i,j] <- log.like.planes(alpha[i],beta[j],y,t)
    }
  }
  like <- exp(loglike - max(loglike))
  main <- expression("Likelihood function for "*alpha*" and "*beta*"")
  contour(alpha,beta,like,drawlabels=F,xlab=expression(alpha),ylab=expression(beta),main=main)
  mle <- expand.grid(alpha,beta)[c(like)==max(c(like)),]
  return(list(alpha=alpha,beta=beta,like=like,mle=mle))
}

newton.raphson.planes <- function(theta0,y,t,abs.tol=1e-6,max.iter=10){
  #################################################################  
  # Run a Newton-Raphson algorithm in the poisson planes model
  #################################################################
  iter <- 1
  theta <- theta0
  f <- loglike.planes(theta,y,t)
  theta.iter <- NULL
  theta.iter <- rbind(theta.iter,t(theta))
  repeat{
    iter <- iter+1
    g <- gradient.planes(theta,y,t)
    H <- hessian.planes(theta,y,t)
    theta.new <- theta - solve(H)%*%g
    f.new <- loglike.planes(theta.new,y,t)
    f.diff <- abs(f.new - f)
    f <- f.new
    theta <- theta.new
    theta.iter <- rbind(theta.iter,t(theta))
    if (f.diff<abs.tol | iter>=max.iter) break
  }
  return(list(theta=theta.iter,mle=theta))
}


# Plot posterior (likelihood) on a grid
par(ask=F)
grid <- contour.likelihood.planes(0,10,0,10,100,accid,time)
grid <- contour.likelihood.planes(20,35,-3,3,100,accid,time)
grid <- contour.likelihood.planes(20,37,-3,1,100,accid,time)
points(grid$mle[1],grid$mle[2],pch="X",col=4,cex=2)
points(lm$coeff[1],lm$coeff[2],col=2,cex=2,pch=16)


# Newton-Raphson for (alpha,beta) using nlminb
theta0 <- c(22,-2)
newton1 <- nlminb(theta0,loglike.planes,y=accid,t=time)
newton2 <- nlminb(theta0,loglike.planes,gradient.planes,hessian.planes,y=accid,t=time)

# Plot the starting values and posterior mode
grid <- contour.likelihood.planes(20,37,-3,1,100,accid,time)
points(theta0[1],theta0[2],pch=16,col=2,cex=2)
points(newton1$par[1],newton1$par[2],pch="X",col=2,cex=2)
points(newton2$par[1],newton2$par[2],pch="X",col=3,cex=2)


# Posterior mode and asymptotic posterior intervals
theta.hat <- newton1$par
covar.hat <- solve(hessian.planes(theta.hat,y=accid,t=time))
alpha.est <- theta.hat[1] + c(-1,0,1)*1.96*sqrt(diag(covar.hat)[1])
beta.est  <- theta.hat[2] + c(-1,0,1)*1.96*sqrt(diag(covar.hat)[2])
alpha.est
beta.est


par(ask=F)
# Show four different starting values 
theta <- rbind(c(22,-2),c(35,-3),c(35,0),c(22,1)) 
grid <- contour.likelihood.planes(20,37,-3,1,100,y=accid,t=time)
points(theta[,1],theta[,2],col=2,pch=16,cex=2)


## Track newton-raphson interations for different starting values
## (function newton.raphson.planes is defined above)
par(ask=T)
for (i in 1:4){
  newton <- newton.raphson.planes(theta[i,],accid,time)
  grid <- contour.likelihood.planes(20,37,-3,1,100,y=accid,t=time)
  points(theta[i,1],theta[i,2],col=2,pch=16,cex=2)
  lines(newton$theta[,1],newton$theta[,2],type="o",col=2,pch=16)
}



#####################################################
########## STAT 289: CODE FOR LECTURE 5a ############
#####################################################

par(mfrow=c(1,1),mar=c(3.5,3.5,2,1),mgp=c(2,.65,0),las=1)

### load library with multivariate normal sampling function
library(MASS)

#############################################
########## READING IN BODYFAT DATA ##########
#############################################

data <- read.table("bodyfat.csv",header=T,sep=",")
attach(data)

#pdf("lecture09-plot1.pdf",height=5,width=8.5)
par(mfrow=c(1,2))
plot(abdomen,weight,pch=19)
plot(wrist,weight,pch=19)
#dev.off()


#############################################
####### MLE FIT OF REGRESSION MODEL #########
#############################################

model <- lm(weight~abdomen+wrist)

#pdf("lecture09-plot2.pdf",height=5,width=8.5)
par(mfrow=c(2,1),mar=c(3.5,3.5,2,1),mgp=c(2,.65,0),las=1)
plot(model$fitted.values,model$residuals,pch=19)
abline(h=0,col=2)
qqnorm(model$residuals)
qqline(model$residuals)
#dev.off()

#############################################
### SAMPLING FROM POSTERIOR DISTRIBUTION ####
#############################################

beta.hat <- model$coef
s <- summary(model)$sigma
s2 <- s^2
V.beta <- summary(model)$cov.unscaled

numsamp <- 1000
numparam <- length(beta.hat)
numobs <- length(wrist)
numdf <- numobs-numparam
beta.samp <- matrix(NA,nrow=numsamp,ncol=numparam)
sigsq.samp <- rep(NA,numsamp)
for (i in 1:numsamp){
  temp <- rchisq(1,numdf)
  cursigsq <- numdf*s2/temp
  curvarbeta <- cursigsq*V.beta
  curbeta <- mvrnorm(1,beta.hat,curvarbeta)
  sigsq.samp[i] <- cursigsq
  beta.samp[i,] <- curbeta
}

#############################################
#### SUMMARIZING POSTERIOR DISTRIBUTIONS ####
#############################################

## posterior means
postmean.beta <- apply(beta.samp,2,mean)
postmean.sigsq <- mean(sigsq.samp)

## posterior correlation between variables
cor(cbind(beta.samp,sigsq.samp))

## 95% posterior intervals
allsamples <- cbind(beta.samp,sigsq.samp)
allsamples.sort <- apply(allsamples,2,sort)
allsamples.sort[c(25,975),]


## posterior histograms
pdf("lecture09-plot3.pdf",height=7,width=7)
par(mfrow=c(2,2),mar=c(3.5,3.5,2,1),mgp=c(2,.65,0),las=1)
hist(allsamples[,1],main="Intercept"); abline(v=c(postmean.beta[1],beta.hat[1]),col=c(2,3))
hist(allsamples[,2],main="Abdomen");   abline(v=c(postmean.beta[2],beta.hat[2]),col=c(2,3))
hist(allsamples[,3],main="Wrist");     abline(v=c(postmean.beta[3],beta.hat[3]),col=c(2,3))
hist(allsamples[,4],main="Sigsq");     abline(v=c(postmean.sigsq,s2),col=c(2,3))
dev.off()


#############################################
#### POSTERIOR PREDICTIVE DISTRIBUTION ######
#### FOR A NEW COVARIATE VECTOR Xstar  ######
#############################################

Xstar <- c(1,130,17)  # new person with abdomen <- 130 and wrist <- 17
Xstar <- t(Xstar)  # making it a row vector


## use posterior samples from before:

ystar.samp <- rep(NA,numsamp)
for (i in 1:numsamp){
   xstarbeta <- Xstar%*%t(t(beta.samp[i,]))
   ystar.samp[i] <- rnorm(1,mean=xstarbeta,sd=sqrt(sigsq.samp[i]))
}   

ystar.postmean <- mean(ystar.samp)

pdf("lecture09-plot4.pdf",height=7,width=6)
par(mfrow=c(2,1))
xmin <- min(weight,ystar.postmean)
xmax <- max(weight,ystar.postmean)
hist(weight,main="Dataset Weights",xlim=c(xmin,xmax))
abline(v=mean(weight),col=2)
hist(ystar.samp,main="Predicted Weight of New Person",xlim=c(xmin,xmax))
abline(v=ystar.postmean,col=2)
dev.off()


#############################################
#### POSTERIOR PREDICTIVE FOR CHECKING ######
####            THE MODEL              ######
#############################################

par(mfrow=c(1,1))
plot(model$fitted.values,model$residuals,pch=19)
plot(model$fitted.values,abs(model$residuals),pch=19)

cor(model$fitted.values,abs(model$residuals))

### take just our first set of sampled values
curbeta <- beta.samp[1,]
cursigsq <- sigsq.samp[1]

### sample whole new set of weights given those parameter values
Xmat <- cbind(rep(1,numobs),abdomen,wrist)

newweight <- rep(NA,numobs)
for (j in 1:numobs){
   Xmatbeta <- Xmat[j,]%*%t(t(curbeta))
   newweight[j] <- rnorm(1,mean=xstarbeta,sd=sqrt(cursigsq))
}   

### for this new vector of weight, re-fit regression model
newmodel <- lm(newweight~abdomen+wrist)
plot(newmodel$fitted.values,newmodel$residuals,pch=19)

Tstat <- cor(newmodel$fitted.values,abs(newmodel$residuals))

### compare to actual data
Tstat.actual <- cor(model$fitted.values,abs(model$residuals))

Tstat 
Tstat.actual

### repeat our posterior predictive sampling for 100 different sampled parameter values
Tstat <- rep(NA,100)

for (i in 1:100){
  newweight <- rep(NA,numobs)
  for (j in 1:numobs){
    Xmatbeta <- Xmat[j,]%*%t(t(beta.samp[j,]))
    newweight[j] <- rnorm(1,mean=xstarbeta,sd=sqrt(sigsq.samp[j]))
  }
  newmodel <- lm(newweight~abdomen+wrist)
  plot(newmodel$fitted.values,newmodel$residuals,pch=19)
  Tstat[i] <- cor(newmodel$fitted.values,abs(newmodel$residuals))
  print(i)
}   

### histogram of posterior predictive dist of Tstat ###
### add line for Tstat for actual data
pdf("lecture09-plot5.pdf",height=5,width=6)
par(mfrow=c(1,1))
hist(Tstat,xlim=range(Tstat,Tstat.actual))
abline(v=Tstat.actual,col=2,lwd=2)
dev.off()

sum(Tstat<Tstat.actual)
