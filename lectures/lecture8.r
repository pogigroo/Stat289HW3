#source("lecture8-aux.r")

#################################################
##### GIBBS FOR NORMAL HIERARCHICAL MODEL #######
#################################################

## simulated data with sigsq=8, tausq=1.5, mu0=0
data <- read.table("normhier.txt")
y<-data
par(mfrow=c(1,1),cex.main=1.5,las=1)
boxplot(y,main="Boxplot of Data by Group") 

true <- c(0,1.5,8)

##Calculating necessary statistics:
m <- ncol(y)
n <- rep(nrow(y),m)
means <- apply(y,2,mean)
ntot <- sum(n)


gibbs1 <- gibbs.norm.hier(y,means,m,n,ntot,1000,c(1,1,1))
gibbs2 <- gibbs.norm.hier(y,means,m,n,ntot,1000,c(10,10,10))

theta1 <- gibbs1$theta; mu1 <- gibbs1$mu
theta2 <- gibbs2$theta; mu2 <- gibbs2$mu

names.theta <- expression(mu[0],tau^2,sigma^2)
names.mu <- expression(mu[1],mu[2],mu[3],mu[4],mu[5],mu[6],mu[7],mu[8],mu[9],mu[10])


### checking convergence 
par(mfrow=c(2,3),cex.main=2,las=1)

for (i in c(1,2,3)){
  plot(theta1[,i],type="l",ylim=range(theta1[,i],theta2[,i]),main=names.theta[i])
  lines(theta2[,i],col=2)
}

for (i in c(1,2,3)){
  plot(mu1[,i],type="l",ylim=range(mu1[,i],mu2[,i]),main=names.mu[i])
  lines(mu2[,i],col=2)
}



### throwing away first 100 as burn-in
good = 101:1000

### autocorrelation of draws in chain 1
par(mfrow=c(2,3),cex.main=2,las=1)
for (i in 1:3) {acf(theta1[good,i],main=names.theta[i])}
for (i in 1:3) {acf(mu1[good,2],main=names.mu[i])}


### thinning chains (every 2nd iteration) and combining chains
good = seq(101,1000,by=2)
theta <- rbind(theta1[good,],theta2[good,])
mu <- rbind(mu1[good,],mu2[good,])


### posterior distributions of parameters (lines for true values)
par(mfrow=c(2,2))
for (i in 1:3){
  hist(theta[,i],xlab="",main=names.theta[i])  
  abline(v=true[i],col=2,lwd=2)
}
apply(theta[,1:3],2,quantile,c(.025,.50,.975))


### Examining Shrinkage
mu.postmean <- apply(mu,2,mean)
mu0.postmean <- mean(theta[,1])
par(mfrow=c(1,1),cex.main=1.5,las=1)
plot(1:10,means,main="Shrinkage of Normal Means",pch=19)
abline(h=mu0.postmean,col=4,lwd=2)
points(1:10,mu.postmean,pch=19,col=2)
legend("topright",c("Data Mean","Mu PostMean","Mu0 PostMean"),pch=19,col=c(1,2,4))


par(mfrow=c(2,1),cex.main=1.5,mar=c(4,4,2,2),las=1)
boxplot(y,main="Original Data (Y)")
boxplot(mu~col(mu),ylim=range(y),main=expression("Posterior Samples of "*mu[i]*""))
round(apply(mu,2,quantile,c(.025,.50,.975)),3)




#################################################
########## GIBBS FOR MIXTURE MODEL ##############
#################################################

#reading in data: alpha=0.25, mu0=0, mu1=3, sigsq0=1, sigsq1=1
data1 <- read.table("normnorm.lec16a.txt")
Y <- data1[,1]
n <- length(Y)
par(mfrow=c(1,1))
hist(Y)

true <- c(.25,0,3,1,1)
names <- expression(alpha,mu[0],mu[1],sigma[0]^2,sigma[1]^2)
gibbs1 <- gibbs.norm.mixture(Y,n,1000,c(.05,0,1,2,2))
params1 <- gibbs1$params


# examining iterations
par(mfrow=c(2,3))
for (i in 1:5){
  plot(params1[,i],type="l",main=names[i])
  abline(h=true[i],col=4)
}


#reading in a less obvious dataset: alpha=0.25, mu0=0, mu1=2, sigsq0=1, sigsq1=1
data2 <- read.table("normnorm.lec16b.txt")
Y <- data2[,1]
n <- length(Y)
par(mfrow=c(1,1))
hist(Y)

true <- c(.25,0,2,1,1)

#running Gibbs sampler for new dataset
gibbs1 <- gibbs.norm.mixture(Y,n,1000,c(.50,0,1,2,2))
gibbs2 <- gibbs.norm.mixture(Y,n,1000,c(.75,1,4,1,1))

params1 <- gibbs1$params
params2 <- gibbs2$params

par(mfrow=c(2,3)) 
for (i in 1:5){
  plot(params1[,i],type="l",main=names[i])
  abline(h=true[i],col=4,lwd=2)
}

par(mfrow=c(2,3)) 
for (i in 1:5){
  plot(params2[,i],type="l",main=names[i])
  abline(h=true[i],col=4,lwd=2)
}


# comparing iterations between two chains 
par(mfrow=c(2,3))
for (i in 1:5){
 plot(params1[,i],type="l",col=2,ylim=range(params1[,i],params2[,i]))
 lines(params2[,i],col=3)
 abline(h=true[i],col=4)
}




# removing 200 iterations of burn-in from each chain
# calculating acf for each chain
good <- 201:1000
par(mfrow=c(2,3)); for (i in 1:5) acf(params1[good,i])
par(mfrow=c(2,3)); for (i in 1:5) acf(params2[good,i])

# need longer chains!  try it with 10000 draws 
# but takes too long for class, so I did it earlier

true <- c(.25,0,2,1,1)
gibbs1 <- gibbs.norm.mixture(Y,n,10000,c(.50,0,1,2,2),plot=T)
gibbs2 <- gibbs.norm.mixture(Y,n,10000,c(.75,1,4,1,1),plot=T)

params1 <- gibbs1$params
params2 <- gibbs2$params

#numiters <- 10000
#params1<-read.table("paramfile.lec16.1.txt")
#params2<-read.table("paramfile.lec16.2.txt")
par(mfrow=c(2,3),cex.main=2)
for (i in 1:5){
  plot(params1[,i],type="l",col=2,ylim=range(params1[,i],params2[,i]),main=names[i])
  lines(params2[,i],col=3)
}


# remove the first 1000 draws as burn-in
good1 <- 1001:10000
par(mfrow=c(2,3)); for (i in 1:5) acf(params1[good1,i],lag.max=100,main=names[i])
par(mfrow=c(2,3)); for (i in 1:5) acf(params2[good1,i],lag.max=100,main=names[i])


# taking only every fiftieth draw 
good2 <- seq(1001,10000,by=50)
par(mfrow=c(2,3)); for (i in 1:5) acf(params1[good2,i],lag.max=100,main=names[i])
par(mfrow=c(2,3)); for (i in 1:5) acf(params2[good2,i],lag.max=100,main=names[i])


# combining chains and calculating posterior intervals
params <- rbind(params1[good1,],params2[good1,])

alpha  <- params[,1]
mu0    <- params[,2]
mu1    <- params[,3]
sigsq0 <- params[,4]
sigsq1 <- params[,5]

par(mfrow=c(1,1))
hist(alpha)
quantile(alpha,c(.025,.5,.975))

post.mean <- apply(params[,1:5],2,mean)
par(mfrow=c(2,3)); for (i in 1:5) {hist(params[,i],main=names[i]); abline(v=post.mean[i],col=2); abline(v=true[i],col=4)}
apply(params[,1:5],2,quantile,c(.025,.5,.975))
post.mean


# fitted densities
numiters <- nrow(params)
par(mfrow=c(1,1),cex.main=1.5,las=1)
hist(Y,prob=T,ylim=c(0,.33))
xx <- seq(-4,6,len=100)
for (i in seq(1,numiters,len=500)){ 
  yy <- alpha[i]*dnorm(xx,mu1[i],sqrt(sigsq1[i]))+(1-alpha[i])*dnorm(xx,mu0[i],sqrt(sigsq0[i]))
  lines(xx,yy,col=2)
}




#################################################
##### NOW APPLY TO BASEBALL DATASET #############
#################################################

#Reading in data (only use data where player-seasons where ab>100)
hitters <- read.table("players.post1970.noP.txt",sep=",",as.is=T,header=T)
hitters <- hitters[hitters[,"AB"]>100,]

#Calculating homerun proportion:
ha <- hitters[,"HR"]/hitters[,"AB"]
hist(ha,main="Histogram of Home Run Average")
n <- length(ha)



#Run the chain from two different starting points
start1 <- c(.5,.01,.1,.01,.01)
start2 <- c(.8,.03,.05,.005,.005)
gibbs1 <- gibbs.norm.mixture(ha,n,2000,start1,plot=F,eval.post=F)
gibbs2 <- gibbs.norm.mixture(ha,n,2000,start2,plot=F,eval.post=F)

params1 <- gibbs1$params;   indics1 <- gibbs1$indics
params2 <- gibbs2$params;   indics2 <- gibbs2$indics

# examining iterations for chains 1 and 2
good <- 100:2000
dev.set(2); par(mfrow=c(2,3),cex.main=2,las=1); for (i in 1:5) plot(params1[good,i],type="l",main=names[i])
dev.set(3); par(mfrow=c(2,3),cex.main=2,las=1); for (i in 1:5) plot(params2[good,i],type="l",main=names[i])


# comparing iterations between two chains 
good1 <- 100:2000
par(mfrow=c(2,3),cex.main=2,las=1)
for (i in 1:5){
  ylim <-  range(params1[good1,i],params2[good1,i]) 
  plot(params2[good1,i],type="l",col=2,ylim=ylim,main=names[i])
  lines(params1[good1,i],col=3)
}

#Label-switching problem, change (alpha,mu1,sigma1) for (1-alpha,mu0,sigma0)
#params0 <- params2[,c(1,3,2,5,4,6)]
#params2 <- params0
#params2[,1] <- 1-params2[,1]


dev.set(2); par(mfrow=c(2,3)); for (i in 1:5) acf(params1[,i],lag.max=100,main=names[i])
dev.set(3); par(mfrow=c(2,3)); for (i in 1:5) acf(params2[,i],lag.max=100,main=names[i])


# taking only every fiftieth draw 
good2 <- seq(200,2000,by=50)
dev.set(2); par(mfrow=c(2,3)); for (i in 1:5) acf(params1[good2,i],main=names[i])
dev.set(3); par(mfrow=c(2,3)); for (i in 1:5) acf(params2[good2,i],main=names[i])


# combining chains and calculating posterior intervals
params <- rbind(params1[good1,],params2[good1,])
indics <- rbind(indics1[good1,],indics2[good1,])
#params <- gibbs2$params[good1,]

par(mfrow=c(2,3))
for (i in 1:5) hist(params[-(1:100),i],xlab="",ylab="",main=names[i],col="skyblue")

dev.set(2)
par(mfrow=c(1,1))
lim <- range(params[,2:3])
plot(params[,2],params[,3],xlim=lim,ylim=lim,xlab=names[2],ylab=names[3],
     main=expression("Joint Posterior for "*mu[0]*" and "*mu[1]*""))
abline(0,1,col=4)

lim <- range(params[,4:5])
plot(params[,4],params[,5],xlim=lim,ylim=lim,xlab=names[4],ylab=names[5],
     main=expression("Joint Posterior for "*sigma[0]^2*" and "*sigma[1]^2*""))
abline(0,1,col=4)


#apply(params[,1:5],2,mean)
#apply(params[,1:5],2,sd)
t(apply(params[,1:5],2,quantile,c(.025,.5,.975)))


alpha  <- params[,1]
mu0    <- params[,2]
mu1    <- params[,3]
sigsq0 <- params[,4]
sigsq1 <- params[,5]

par(mfrow=c(1,1))
hist(alpha)
quantile(alpha,c(.025,.5,.975))

#post.mean <- apply(params[,1:5],2,mean)
#par(mfrow=c(2,3)); for (i in 1:5) {hist(params[,i],main=names[i]); abline(v=post.mean[i],col=2); abline(v=true[i],col=4)}
#apply(params[,1:5],2,quantile,c(.025,.5,.975))
#post.mean



## fitted densities
numiters <- nrow(params)
par(mfrow=c(1,1),cex.main=1.5)
hist(ha,prob=T,ylim=c(0,27),main="Histogram of Home Run Averages")
xx <- seq(0,.16,len=100)
for (i in seq(1,numiters,len=500)){ 
  yy <- alpha[i]*dnorm(xx,mu1[i],sqrt(sigsq1[i]))+(1-alpha[i])*dnorm(xx,mu0[i],sqrt(sigsq0[i]))
  lines(xx,yy,col=3)
}

## identify Bobby Abreu's seasons
abreu <- which(hitters[,"name"]=="abreubo01")

## posterior probabilities that Abreu is elite
gamma <- apply(indics[,abreu],2,mean)
cbind(hitters[abreu,c("name","year","AB","HR")],ha[abreu],gamma)

abline(v=ha[abreu],col=2:9,lwd=1)
legend(.10,25,lty=1,col=2:9,lwd=2,paste("Abreu",1997:2004))#,"(",round(gamma,3),")")
