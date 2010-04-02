############################################################
############## R CODE FOR STAT 289 LECTURE 4 ###############
############################################################
 
par(mfrow=c(1,1),mar=c(3.5,3.5,2,1),mgp=c(2,.65,0),las=1)

summarize.discrete = function(x,p){
  ###########################################################
  # Auxiliary function:
  # returns mode, mean, sd of discrete distribution p(x)
  ###########################################################
  p = p/sum(p)
  mode = x[which(p==max(p))]
  mean = sum(x*p)
  sdev = sqrt(sum(x^2*p)-mean^2)
  return(list(mode=mode,mean=mean,sd=sdev))
}



###################################################################
#### 1. Baseball Data
#### Grid Search and Grid Sample: Semi-Conjugate Normal Example ###
###################################################################

## Reading in Data:
data <- read.table("players.post1970.txt",header=T,sep=",")

## Remove pitchers, player-seasons where ab<=100, and NA rows.
subset <- data$pos!="P" & data$AB>100 & !is.na(data[,1])
data <- data[subset,]

## Calculating batting average
ba <- data$H/data$AB
n <- length(ba)

hist(ba,xlab="Batting Average (y)",main="MLB Batting Averages, 1970-2004")

min(ba)
data[which(ba==min(ba)),]
max(ba)
data[which(ba==max(ba)),]




eval.post.sigsq = function(sigsq,y,mu0,nu0,tausq0,sigsq0){
  ##########################################################
  ## Marginal posterior for sigma^2 in semi-conjugate model
  ##########################################################
  n = length(y)
  ybar = mean(y)
  postvar = 1/(1/tausq0 + n/sigsq0)
  postmean = postvar*(mu0/tausq0 + ybar*n/sigsq0)
  deviance = sum((y-postmean)^2)
  
  logpost = 0.5*log(postvar) + dnorm(postmean,mu0,sqrt(tausq0),log=T) 
  logpost = logpost -((nu0+n)/2+1)*log(sigsq) - (nu0*sigsq0+deviance)/(2*sigsq)
  post = exp(logpost-max(logpost))
  return(post)
}




## setting hyperparameter values
tausq0 <- 0.001 
sigsq0 <- 0.001 
nu0 <- 10
mu0 <- 0.2


## Evaluate the marginal posterior for sigmasq

type = "l"

sigmasq <- seq(1e-6,5e-3,len=100)
post.sigmasq <- eval.post.sigsq(sigmasq,ba,mu0,nu0,tausq0,sigsq0)
plot(sigmasq,post.sigmasq,type=type,xlab="sigma^2",main="Posterior Distribution of sigma^2")

sigmasq <- seq(.001,.002,len=100)
post.sigmasq <- eval.post.sigsq(sigmasq,ba,mu0,nu0,tausq0,sigsq0)
plot(sigmasq,post.sigmasq,type=type,xlab="sigma^2",main="Posterior Distribution of sigma^2")


sigmasq <- seq(.00115,.0013,len=100)
post.sigmasq <- eval.post.sigsq(sigmasq,ba,mu0,nu0,tausq0,sigsq0)
post.sigmasq <- post.sigmasq/sum(post.sigmasq)
a <- summarize.discrete(sigmasq,post.sigmasq)

plot(sigmasq,post.sigmasq,type="h",xlab="sigma^2",main="Posterior Dist. of sigma^2 (Semi-Conjugate Prior)")
abline(v=var(ba),lty=2,col=4)
abline(v=a$mode,lty=1,col=2)
legend("topright",lty=c(2,1),col=c(4,2),c("sample variance","posterior mode"))






## grid sampling: sample 1000 values of sigmasq proportional to post.sigmasq

sigmasq.samp <- sample(sigmasq,size=1000,replace=T,prob=post.sigmasq)

plot(sigmasq,post.sigmasq,type="l",xlab="sigma^2",xlim=range(sigmasq),ylab="density",
     col=2,lwd=2,main="Posterior Samples of sigma^2")
par(new=T)
hist(sigmasq.samp,prob=T,xlab="sigma^2",ylab="",xlim=range(sigmasq))
#for (i in 1:5) text(.0013,26000-1500*i,txt[i],adj=1)
quantile(sigmasq.samp,c(.025,.5,.975))




## sampling mu given sampled sigmasq

postvar <- 1/(1/tausq0 + n/sigmasq.samp)
postmean <- postvar*(mu0/tausq0 + mean(ba)*n/sigmasq.samp)
mu.samp <- rnorm(1000,postmean,sqrt(postvar))

hist(mu.samp,xlab="mu",main="Posterior Samples of mu (Semi-Conjugate Prior)")
abline(v=mean(ba),lty=2,col=4)
abline(v=mean(mu.samp),lty=1,col=2)
legend("topleft",lty=c(2,1),col=c(4,2),c("sample mean","posterior mean"))





################################################################
### 2. Grid Search and Grid Sample: Poisson Planes Dataset #####
################################################################

## input data:  
data <- read.table("planes.txt",skip=1)
y <- data[,2]
t <- data[,1]-1976
n <- length(y)


par(mfrow=c(2,1),mar=c(4,4,2,1),las=1)
hist(y,xlim=c(10,40),main="Annual Flight Fatalities, 1976-1985"); box()
plot(t,y,main="Fatalities vs Year",pch=16)


log.like.planes = function(alpha,beta,y,t){
  ###########################################################
  ## Loglikelihood function for modified poisson regression
  ###########################################################
  if (alpha+beta*max(t)<=0) value = -Inf
  if (alpha+beta*max(t)>0) value = sum(y*log(alpha+beta*t)-(alpha+beta*t))
  return(value)
}



# Evaluate the posterior (=likelihood) distribution over a grid
ngrid = 100
if (0){
 alpha <- seq(25,35,len=ngrid)
 beta <- seq(-2,0,length=ngrid)
}
alpha <- seq(20,37,len=ngrid)
beta <- seq(-3,1,length=ngrid)

logpost <- matrix(NA,ngrid,ngrid)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    logpost[i,j] <- log.like.planes(alpha[i],beta[j],y,t)
  }
}

post <- exp(logpost-max(logpost))
post <- post/sum(post)
sum(post)


# Contour plot of the joint posterior
par(mfrow=c(1,1))
contour(alpha,beta,post,xlab="alpha",ylab="beta",drawlabels=F,main="Joint Posterior of alpha and beta",col=2,lwd=2)
#points(expand.grid(alpha,beta),pch=16,cex=.5,col="gray60")



## Computing marginal posterior of alpha
par(mfrow=c(2,1))
post.alpha <- apply(post,1,sum)
plot(alpha,post.alpha,type="l",ylab="",main="Marginal Posterior of alpha",col=4,lwd=2)
a <- summarize.discrete(alpha,post.alpha)
a


# Computing marginal posterior of beta
post.beta <- apply(post,2,sum)
plot(beta,post.beta,type="l",ylab="",main="Marginal Posterior of beta",col=4,lwd=2)
b <- summarize.discrete(beta,post.beta)
b


# Computing conditional posterior of beta for different values of alpha
post.beta.cond <- matrix(NA,ngrid,ngrid)
for (i in 1:ngrid) post.beta.cond[i,] <- post[i,]/post.alpha[i]


# Plot conditional posterior of beta for different alphas
index = c(.05,.5,.95)*ngrid
alpha[index]
par(mfrow=c(3,1))
for (i in index){
  title <- paste("Conditional Posterior of beta for alpha=",round(alpha[i],2))
  plot(beta,post.beta.cond[i,],type="l",ylab="",main=title,col=4,lwd=2)
}



## Sampling grid values 
alpha.samp <- rep(NA,10000)
beta.samp <- rep(NA,10000)
for (m in 1:10000){
  a <- sample(1:ngrid,size=1,prob=post.alpha)
  b <- sample(1:ngrid,size=1,prob=post.beta.cond[a,])
  alpha.samp[m] <- alpha[a]
  beta.samp[m] <- beta[b]
}


par(mfrow=c(1,1))
contour(alpha,beta,post,drawlabels=F,col=2,main="Posterior Samples of alpha and beta")
points(alpha.samp,beta.samp,pch=16,cex=.35)
points(jitter(alpha.samp),jitter(beta.samp),pch=16,cex=.35)  # Jitter points to show density better



## calculating posterior means/intervals for alpha and beta
par(mfrow=c(2,1))
hist(alpha.samp,main="Posterior Samples of alpha")
hist(beta.samp,main="Posterior Samples of beta")


mean(alpha.samp)
mean(beta.samp)
quantile(alpha.samp,c(.025,.975))
quantile(beta.samp,c(.025,.975))
mean(beta.samp>0)



## Plot posterior samples of alpha+beta*x
plot(t,y,main="Posterior Samples of alpha+beta*x")
for (i in 1:1000) abline(alpha.samp[i],beta.samp[i],col=3)
points(t,y,pch=16)


## Predicted new observation for 1986 (t = 10):
pred.rate <- alpha.samp + beta.samp*10
pred.accidents <- rpois(10000,pred.rate)
hist(pred.accidents,main="Posterior Predictive Distribution for Accidents at t=10")

mean(pred.accidents)
quantile(pred.accidents,c(.025,.975))
mean(pred.accidents>30)




############################################################
### 3. Grid Search and Grid Sample: Bioassay Dataset #######
############################################################

## input data:  
data <- read.table("bioassay.txt",header=T)
x <- data[,1]
n <- data[,2]
y <- data[,3]
N <- length(y)

plot(x,y/n,xlab="Dose",ylab="prop.fatal",pch=19,main="Proportion of Fatalities vs Dose")


## logistic regression estimates of alpha and beta:
mat <- cbind(y,n-y)
glm <- glm(mat~x,family="binomial")
summary(glm)





log.like.logistic <- function(alpha,beta,y,x,n,N){
  #####################################################
  ## Loglikelihood function for logistic regression
  #####################################################
  theta <- exp(alpha+beta*x)/(1+exp(alpha+beta*x))
  loglike <- sum(y*log(theta) + (n-y)*log(1-theta))
  return(loglike)
}



## Evaluate joint posterior over a range of alpha and beta
ngrid <- 100
alpha <- seq(0,1,len=ngrid)
beta <- seq(0,10,len=ngrid)
logpost <- matrix(NA,ngrid,ngrid)

for (i in 1:ngrid){
  for (j in 1:ngrid){
    logpost[i,j] <- log.like.logistic(alpha[i],beta[j],y,x,n,N)
  }
}
post <- exp(logpost-max(logpost))
post <- post/sum(post)

contour(alpha,beta,post,drawlabels=F,main="Joint Posterior of alpha and beta")  



## Evaluate joint posterior over a new range of alpha and beta
ngrid <- 100
alpha <- seq(-2,4,len=ngrid)
beta <- seq(0,30,len=ngrid)
logpost <- matrix(NA,ngrid,ngrid)

for (i in 1:ngrid){
  for (j in 1:ngrid){
    logpost[i,j] <- log.like.logistic(alpha[i],beta[j],y,x,n,N)
  }
}
post <- exp(logpost-max(logpost))
post <- post/sum(post)


## Contour plot of the joint posterior of alpha and beta
par(mfrow=c(1,1))
contour(alpha,beta,post,drawlabels=F,main="Joint Posterior of alpha and beta")
points(glm$coeff[1],glm$coeff[2],pch="X",col=2,cex=1.5)



par(mfrow=c(2,1))
post.alpha <- apply(post,1,sum)
plot(alpha,post.alpha,type="l",ylab="",main="Marginal Posterior of alpha",col=4,lwd=2)
a <- summarize.discrete(alpha,post.alpha)
a


# Marginal posterior of beta
post.beta <- apply(post,2,sum)
plot(beta,post.beta,type="l",ylab="",main="Marginal Posterior of beta",col=4,lwd=2)
b <- summarize.discrete(beta,post.beta)
b




# Conditional posterior of beta given alpha
post.beta.cond <- matrix(NA,ngrid,ngrid)
for (i in 1:ngrid) post.beta.cond[i,] <- post[i,]/post.alpha[i]

index = c(.20,.80)*ngrid
alpha[index]
par(mfrow=c(2,1))
for (i in index){
  title <- paste("Conditional Posterior of beta for alpha=",round(alpha[i],2))
  plot(beta,post.beta.cond[i,],type="l",ylab="",main=title,col=4,lwd=2)
}





## Sampling grid values from the joint distribution
alpha.beta <- as.matrix(expand.grid(alpha,beta))
post.wts <- c(post)
index <- sample(1:ngrid^2,10000,replace=T,prob=post.wts)
alpha.samp <- alpha.beta[index,1]
beta.samp <- alpha.beta[index,2]

par(mfrow=c(3,1))
contour(alpha,beta,post,drawlabels=F,col=2,main="Posterior Samples of alpha and beta")
points(alpha.samp,beta.samp,pch=16,cex=.35)
points(jitter(alpha.samp),jitter(beta.samp),pch=16,cex=.35)


## calculating posterior means/intervals for alpha and beta
par(mfrow=c(2,1))
hist(alpha.samp,prob=T,main="Posterior Samples of alpha")
text(3,.350,paste("E(alpha|y)=",round(mean(alpha.samp),3)))
text(3,.325,paste("SD(alpha|y)=",round(sd(alpha.samp),4)))

hist(beta.samp,prob=T,main="Posterior Samples of beta")
text(24,.070,paste("E(beta|y)=",round(mean(beta.samp),3)))
text(24,.065,paste("SD(beta|y)=",round(sd(beta.samp),4)))

mean(alpha.samp)
mean(beta.samp)
sd(alpha.samp)
sd(beta.samp)
quantile(alpha.samp,c(.025,.975))
quantile(beta.samp,c(.025,.975))
mean(beta.samp<0)



## Draw posterior samples of theta = exp(a+bx)/(1+exp(a+bx))
par(mfrow=c(1,1))
plot(x,y/n,main="Posterior Samples of Theta(x)")
for (i in 1:1000){
  xgrid <- seq(-1,1,len=100)
  temp <- exp(alpha.samp[i]+beta.samp[i]*xgrid)
  prob <- temp/(1+temp)
  lines(xgrid,prob,col=3)
}
points(x,y/n,pch=16)





## Predicted number of deaths for a new batch with n=5 animals and dosage x=0
par(mfrow=c(2,1))
xnew = 0
pred.rate <- exp(alpha.samp + beta.samp*xnew)/(1+exp(alpha.samp + beta.samp*xnew))
pred.deaths <- rbinom(10000,5,pred.rate)
hist(pred.deaths,xlim=c(0,5),ylab="",main=paste("Posterior Predictive for Deaths when n=5, x=",xnew))

mean(pred.deaths)
quantile(pred.deaths,c(.025,.975))
mean(pred.deaths==5)



## Predicted number of deaths for a new batch with n=5 animals and dosage x=0.5 

xnew = 0.5
pred.rate <- exp(alpha.samp + beta.samp*xnew)/(1+exp(alpha.samp + beta.samp*xnew))
pred.deaths <- rbinom(10000,5,pred.rate)
hist(pred.deaths,xlim=c(0,5),ylab="",main=paste("Posterior Predictive for Deaths when n=5, x=",xnew))

mean(pred.deaths)
quantile(pred.deaths,c(.025,.975))
mean(pred.deaths==5)






#####################################################
############### STAT 289, LECTURE 4A ################
#####################################################

#############################################
## UNIVARIATE NEWTON-RAPHSON: GAMMA DATA ####
#############################################

## input data:  (gamma data with true alpha = 2 and beta = 3)
y <- read.table("gammaexample.txt")
y <- y[,2]
n <- length(y)


## graphing likelihood over range of alpha and beta
loglik.gamma <- function(alpha,beta,x){
  n <- length(x)
  sumx <- sum(x)
  sumlogx <- sum(log(x))
  loglik <- n*alpha*log(beta)-n*lgamma(alpha)
  loglik <- loglik + (alpha-1)*sumlogx - beta*sumx
  loglik
}

alpharange <- ppoints(100)*3+1.5   # alpha between 1.5 and 4.5
betarange <- ppoints(100)*4+2  # beta between 2 and 6
z <- matrix(NA,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    z[i,j] <- loglik.gamma(alpharange[i],betarange[j],y)
  }
}
z <- exp(z-max(z))
contour(alpharange,betarange,z,xlab="alpha",ylab="beta",drawlabels=F)



## Plotting the "Score" Function:
alpharange <- ppoints(100)*3+1.5   # alpha between 1.5 and 4.5
frange <- n*log(3)+sum(log(y))-n*digamma(alpharange)
plot(alpharange,frange,type="l")
abline(h=0)



## Newton-Raphson algorithm
gamma.1D.NR <- function(old){
  new <- old + (n*log(3)+sum(log(y))-n*digamma(old))/(n*trigamma(old))
  new
}   

iters1 <- c(NA,20); iters1[1] <- 3.00  # starting value
iters2 <- c(NA,20); iters2[1] <- 0.01  # starting value
iters3 <- c(NA,20); iters3[1] <- 10.0  # starting value

for (i in 1:19) {iters1[i+1] <- gamma.1D.NR(iters1[i])}
for (i in 1:19) {iters2[i+1] <- gamma.1D.NR(iters2[i])}
for (i in 1:19) {iters3[i+1] <- gamma.1D.NR(iters3[i])}


plot(iters1,type="n",ylim=range(iters1,iters2,iters3))
lines(iters1,type="o",col=1,pch=15,cex=.75)
lines(iters2,type="o",col=2,pch=16,cex=.75)
lines(iters3,type="o",col=3,pch=17,cex=.75)



#############################################
## MULTIVARIATE NEWTON-RAPHSON: GAMMA DATA ##
#############################################

## input data:  (gamma data with true alpha = 2 and beta = 3)
y <- read.table("gammaexample.txt")
y <- y[,2]
n <- length(y)

## graphing likelihood over range of alpha and beta
loglik.gamma <- function(alpha,beta,x){
  n <- length(x)
  sumx <- sum(x)
  sumlogx <- sum(log(x))
  loglik <- n*alpha*log(beta)-n*lgamma(alpha)
  loglik <- loglik + (alpha-1)*sumlogx - beta*sumx
  loglik
}
alpharange <- ppoints(100)*3+1.5   # alpha between 1.5 and 4.5
betarange <- ppoints(100)*4+2  # beta between 2 and 6
z <- matrix(NA,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    z[i,j] <- loglik.gamma(alpharange[i],betarange[j],y)
  }
}
z <- exp(z-max(z))
contour(alpharange,betarange,z,xlab="alpha",ylab="beta",drawlabels=F)


## Multivariate Newton Raphson algorithm:

gamma.2D.NR <- function(olda,oldb,y){
  n <- length(y)
  Jaa <- -n*trigamma(olda)
  Jab <- n/oldb
  Jbb <- -n*olda/(oldb^2)
  J <- rbind(c(Jaa,Jab),c(Jab,Jbb))
  Jinv <- solve(J)
  ga <- n*log(oldb) - n*digamma(olda) + sum(log(y))
  gb <- n*olda/oldb - sum(y)
  gvec <- t(t(c(ga,gb)))
  oldvec <- t(t(c(olda,oldb)))
  newvec <- oldvec - Jinv%*%gvec
  newvec
}
iters <- matrix(NA,nrow=20,ncol=2)
iters[1,] <- c(2,3)  # starting value
for (i in 1:19){
  new <- gamma.2D.NR(iters[i,1],iters[i,2],y)
  iters[i+1,] <- new
}

contour(alpharange,betarange,z,xlab="alpha",ylab="beta",drawlabels=F)
lines(iters,type="o",col=4,pch=16)


## better way to do iteration:
#old <- c(2,3) # starting value
old <- c(2.5,2.5) # starting value
itersmat <- NULL
itersmat <- rbind(itersmat,old)
new <- c(0,0)
diff <- sum(abs(new-old))
while (diff > 0.0000001){
  new <- gamma.2D.NR(old[1],old[2],y)
  itersmat <- rbind(itersmat,c(new[1],new[2]))
  diff <- sum(abs(new-old))
  old <- new
}  


contour(alpharange,betarange,z,xlab="alpha",ylab="beta",drawlabels=F)
lines(itersmat,type="o",col=2,pch=16)



#############################################
# MULTIVARIATE NEWTON-RAPHSON: POISSON DATA #
#############################################

## input data:  
data <- read.table("planes.txt",skip=1)
accid <- data[,2]
time <- data[,1]-1976
n <- length(y)
plot(time,accid)


## graphing likelihood over range of alpha and beta:

loglik.pois <- function(alpha,beta,data.y,data.t){
  n <- length(data.y)
  loglik <- sum(data.y*log(alpha+beta*data.t))
  loglik <- loglik - n*alpha - beta*sum(data.t)
  loglik
}

alpharange <- ppoints(100)*10   # alpha between 0 and 10
betarange <- ppoints(100)*10  # beta between 0 and 10
z <- matrix(-Inf,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    if (alpharange[i] + betarange[j]*max(time) > 0){  # makes sure lambda is positive
      z[i,j] <- loglik.pois(alpharange[i],betarange[j],accid,time)
    }
  }
}
z <- exp(z - max(z))
contour(alpharange,betarange,z,xlab="alpha",ylab="beta",drawlabels=F)



alpharange <- ppoints(100)*15+20   # alpha between 20 and 35
betarange <- ppoints(100)*6-3  # beta between -3 and 3
z <- matrix(-Inf,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    if (alpharange[i] + betarange[j]*max(time) > 0){ # makes sure lambda is positive
      z[i,j] <- loglik.pois(alpharange[i],betarange[j],accid,time)
    }
  }
}
z <- exp(z - max(z))
contour(alpharange,betarange,z,xlab="alpha",ylab="beta",drawlabels=F)


## Multivariate Newton Raphson algorithm 

planes.NR <- function(olda,oldb,data.y,data.t){
  x <- olda+oldb*data.t
  x2 <- x^2
  n <- length(data.y)
  Jaa <- -sum(data.y/x2)
  Jab <- -sum(data.y*data.t/x2)
  Jbb <- -sum(data.y*(data.t^2)/x2)
  J <- rbind(c(Jaa,Jab),c(Jab,Jbb))
  Jinv <- solve(J,tol=1e-14)
  ga <- sum(data.y/x) - n
  gb <- sum((data.y*data.t)/x) - sum(data.t)
  gvec <- t(t(c(ga,gb)))
  oldvec <- t(t(c(olda,oldb)))
  newvec <- oldvec - Jinv%*%gvec
  newvec
}

old <- c(25,-1) # starting value
itersmat1 <- NULL
itersmat1 <- rbind(itersmat1,old)
new <- c(0,0)
diff <- sum(abs(new-old))
while (diff > 0.0000001){
  new <- planes.NR(old[1],old[2],accid,time)
  itersmat <- rbind(itersmat1,c(new[1],new[2]))
  diff <- sum(abs(new-old))
  old <- new
}  


old <- c(20,1) # different starting value
itersmat2 <- NULL
itersmat2 <- rbind(itersmat2,old)
new <- c(0,0)
diff <- sum(abs(new-old))
while (diff > 0.0000001){
  new <- planes.NR(old[1],old[2],accid,time)
  itersmat2 <- rbind(itersmat2,c(new[1],new[2]))
  diff <- sum(abs(new-old))
  old <- new
}  

#old <- c(100,-3) # different starting value
old <- c(30,-2) # different starting value
itersmat3 <- NULL
itersmat3 <- rbind(itersmat3,old)
new <- c(0,0)
diff <- sum(abs(new-old))
while (diff > 0.0000001){
  new <- planes.NR(old[1],old[2],accid,time)
  itersmat3 <- rbind(itersmat3,c(new[1],new[2]))
  diff <- sum(abs(new-old))
  old <- new
}  




contour(alpharange,betarange,z,xlab="alpha",ylab="beta",drawlabels=F)
lines(itersmat1,type="o",col=2,pch=16)
lines(itersmat2,type="o",col=3,pch=16)
lines(itersmat3,type="o",col=4,pch=16)


