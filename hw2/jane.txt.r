12
data<-read.table("planes.txt",skip=1)
y<-data[,2]
t<-data[,1]
n<-length(y)

## (b)
sigma1<-9
sigma2<-0.16
mu1<-28
mu2<-0
pho=-.5
gridx<-seq(20,35,length=1000)
gridy<-seq(-1,1,length=1000)
likelihood<-function(gridx,gridy){
value<-((gridx-mu1)^2/sigma1+(gridy-mu2)^2/sigma2-2*pho*(gridx-mu1)*(gridy-mu2)/(sigma1+sigma2))/(-2*(1-pho^2))
return(value)}
likeli<-matrix(NA,1000,1000)
for (i in 1:1000){
  for (j in 1:1000){
    likeli[i,j]<-likelihood(gridx[i],gridy[j])
  }
}
like<-exp(likeli-exp(likeli))
like<-like/sum(like)
sum(like)
contour(gridx,gridy,like,drawlabels=F)

(e)
result<-lm(y~t)
plot(t,y,main="Linear Regression of alpha and beta")
abline(result)
summary(result)

t1<-t-1975
result1<-lm(y~t1)
summary(result1)

(f)
data<-read.table("planes.txt",skip=1)
y<-data[,2]
t<-data[,1]
t1<-t-1975
n<-length(y)
log.like<-function(y,t1,alpha,beta){
    if (beta*max(t1)+alpha<=0) value=-Inf
    if (beta*max(t1)+alpha>0) value<-sum(y*log(beta*t1+alpha)-(beta*t1+alpha))
    return(value)
}
ngrid<-1000
alpha<-seq(20,40,length=ngrid)
beta<-seq(-3,1,length=ngrid)
loglike<-matrix(NA,ngrid,ngrid)
for (i in 1:ngrid) {
   for (j in 1:ngrid) {
      loglike[i,j]<-log.like(y,t1,alpha[i],beta[j])
   }
}
like<-exp(loglike-max(loglike))
like<-like/sum(like)
sum(like)
contour(alpha,beta,like,xlab="alpha",ylab="beta",drawlabels=F,main="contour of the joint posterior density of (alpha,beta)",col=4)

(g)
post.alpha<-apply(like,1,sum)
post.beta.cond<-matrix(NA,ngrid,ngrid)
for (i in 1:ngrid) post.beta.cond[i,]<-like[i,]/post.alpha[i]
alpha.samp<-rep(NA,ngrid)
beta.samp<-rep(NA,ngrid)
for (m in 1:ngrid){
  a<-sample(1:ngrid,size=1,prob=post.alpha)
  b<-sample(1:ngrid,size=1,prob=post.beta.cond[a,])
  alpha.samp[m]<-alpha[a]
  beta.samp[m]<-beta[b]
}

(h)
pred.rate<-alpha.samp+beta.samp*11
pred.accidents<-rpois(1000,pred.rate)
hist(pred.accidents,main="Histogram of Predictive accidents in 1986")
lower=mean(pred.accidents)-1.96*sqrt(var(pred.accidents))
lower
upper=mean(pred.accidents)+1.96*sqrt(var(pred.accidents))
upper


Question 3
data <- read.table("hitters.post1975.txt",header=T,sep=",")
subset<- data$AB>100 & data$HR>1
data<-data[subset,]
h<-data[,10]
hr<-data[,13]
ab<-data[,8]
n<-length(h)
data<-matrix(c(h,hr,ab),nrow=length(h),byrow=FALSE)
bat_average<-rep(0,n)
for (i in 1:n){
if (ab[i]>0){bat_average[i]<-h[i]/ab[i]}}
par(mfrow=c(1,1))
hist(bat_average,main="Histogram of batting average")
mean.bat.average<-mean(bat_average)
max(bat_average)
which.max(bat_average)
min(bat_average)
which.min(bat_average)
mean(bat_average)

Question 4
log.like<-function(y,n,a,b){
   value=n*lgamma(a+b)-n*lgamma(a)-n*lgamma(b)+sum((a-1)*log(y)+(b-1)*log(1-y))
   return(value)
}
n<-length(bat_average)
ngrid<-100
a<-seq(45,60,len=ngrid)
b<-seq(120,170,len=ngrid)
loglike<-matrix(NA,ngrid,ngrid)
for (i in 1:ngrid){
   for(j in 1:ngrid){
      loglike[i,j]<-log.like(bat_average,n,a[i],b[j])
   }
}
contour(a,b,loglike,xlab="a",ylab="b",drawlabels=F,main="Contour of loglikelihood")
mle<-expand.grid(a,b)[c(loglike)==max(c(loglike)),]

## contour plot using 100 random points from bat_average
log.like<-function(y,n,a,b){
   value=n*lgamma(a+b)-n*lgamma(a)-n*lgamma(b)+sum((a-1)*log(y)+(b-1)*log(1-y))
   return(value)
}
n<-length(bat_average)
ngrid<-100
a<-seq(45,60,len=ngrid)
b<-seq(120,170,len=ngrid)
loglike<-matrix(NA,ngrid,ngrid)
for (i in 1:ngrid){
   for(j in 1:ngrid){
      loglike[i,j]<-log.like(sample(bat_average,100),100,a[i],b[j])
   }
}
contour(a,b,loglike,xlab="a",ylab="b",drawlabels=F,main="Contour of loglikelihood")
mle<-expand.grid(a,b)[c(loglike)==max(c(loglike)),]


like<-exp(loglike-max(loglike))
like<-like/sum(like)
sum(like)

##Question 5
mle<-expand.grid(a,b)[c(loglike)==max(c(loglike)),]
## use 100 random points from bat_average
log.like<-function(y,n,a,b){
   value=n*lgamma(a+b)-n*lgamma(a)-n*lgamma(b)+sum((a-1)*log(y)+(b-1)*log(1-y))
   return(value)
}
n<-length(bat_average)
ngrid<-100
a<-seq(45,60,len=ngrid)
b<-seq(120,170,len=ngrid)
loglike<-matrix(NA,ngrid,ngrid)
for (i in 1:ngrid){
   for(j in 1:ngrid){
      loglike[i,j]<-log.like(sample(bat_average,100),100,a[i],b[j])
   }
}
contour(a,b,loglike,xlab="a",ylab="b",drawlabels=F,main="Contour of loglikelihood")
mle<-expand.grid(a,b)[c(loglike)==max(c(loglike)),]

## Question 6
log.like<-function(y,n,theta){
   value=n*lgamma(theta[1]+theta[2])-n*lgamma(theta[1])-n*lgamma(theta[2])+sum((theta[1]-1)*log(y)+(theta[2]-1)*log(1-y))
   return(value)
}
gradient<-function(theta,n,y){
g=rep(NA,2)
g[1]=n*(digamma(theta[1])-digamma(theta[1]+theta[2]))+sum(log(y))
g[2]=n*(digamma(theta[2])-digamma(theta[1]+theta[2]))+sum(log(1-y))
return(-g)
}
hessian<-function(theta,n,y){
h=matrix(NA,2,2)
h[1,1]=n*(trigamma(theta[1])-trigamma(theta[1]+theta[2]))
h[1,2]=-n*(trigamma(theta[1]+theta[2]))
h[2,1]=-n*(trigamma(theta[1]+theta[2]))
h[2,2]=n*(trigamma(theta[2])-trigamma(theta[1]+theta[2]))
return(-h)
}
####newton=function(theta0,y,abs.tol=1e-6,max.iter=2){
iter=1
theta=theta0
f=log.like(theta,n,y)
repeat{
   iter=iter+1
   g=gradient(theta,n,y)
   h=hessian(theta,n,y)
   theta.new=theta-solve(h)%*%g
   f.new=log.like(theta.new,n,y)
   f.diff=abs(f.new-f)
   f=f.new
   theta=theta.new
   if(f.diff<abs.tol|iter>=max.iter) break
}
return(list(mle=theta))
}####
theta0<-c(50,150)
newton=nlminb(theta0,log.like,gradient,hessian,y=bat_average,n=length(bat_average))

## Question 7

##histogram of home run average
data <- read.table("hitters.post1975.txt",header=T,sep=",")
subset<- data$AB>100 & data$HR>1
data<-data[subset,]
h<-data[,10]
hr<-data[,13]
ab<-data[,8]
n<-length(h)
data<-matrix(c(h,hr,ab),nrow=length(h),byrow=FALSE)
hr_average<-rep(0,n)
for (i in 1:n){
if (ab[i]>0){hr_average[i]<-hr[i]/ab[i]}}
par(mfrow=c(1,1))
hist(hr_average,main="Histogram of home run average")
mean.hr.average<-mean(hr_average)
max(hr_average)
which.max(hr_average)
min(hr_average)
which.min(hr_average)
mean(hr_average)

## contour of loglikelohood of home run average
log.like<-function(y,n,a,b){
   value=n*lgamma(a+b)-n*lgamma(a)-n*lgamma(b)+sum((a-1)*log(y)+(b-1)*log(1-y))
   return(value)
}
n<-length(hr_average)
ngrid<-100
a<-seq(1,7,len=ngrid)
b<-seq(60,150,len=ngrid)
loglike<-matrix(NA,ngrid,ngrid)
for (i in 1:ngrid){
   for(j in 1:ngrid){
      loglike[i,j]<-log.like(hr_average,n,a[i],b[j])
   }
}
contour(a,b,loglike,xlab="a",ylab="b",drawlabels=F,main="Contour of loglikelihood")

## MLE of a and b
mle<-expand.grid(a,b)[c(loglike)==max(c(loglike)),]

## contour plot using 100 random points from hr_average
log.like<-function(y,n,a,b){
   value=n*lgamma(a+b)-n*lgamma(a)-n*lgamma(b)+sum((a-1)*log(y)+(b-1)*log(1-y))
   return(value)
}
n<-length(hr_average)
ngrid<-100
a<-seq(1,7,len=ngrid)
b<-seq(60,150,len=ngrid)
loglike<-matrix(NA,ngrid,ngrid)
for (i in 1:ngrid){
   for(j in 1:ngrid){
      loglike[i,j]<-log.like(sample(hr_average,100),100,a[i],b[j])
   }
}
contour(a,b,loglike,xlab="a",ylab="b",drawlabels=F,main="Contour of loglikelihood")
mle<-expand.grid(a,b)[c(loglike)==max(c(loglike)),]

## Newton-Raphson MLE of a and b
log.like<-function(y,n,theta){
   value=n*lgamma(theta[1]+theta[2])-n*lgamma(theta[1])-n*lgamma(theta[2])+sum((theta[1]-1)*log(y)+(theta[2]-1)*log(1-y))
   return(value)
}
gradient<-function(theta,n,y){
g=rep(NA,2)
g[1]=n*(digamma(theta[1])-digamma(theta[1]+theta[2]))+sum(log(y))
g[2]=n*(digamma(theta[2])-digamma(theta[1]+theta[2]))+sum(log(1-y))
return(-g)
}
hessian<-function(theta,n,y){
h=matrix(NA,2,2)
h[1,1]=n*(trigamma(theta[1])-trigamma(theta[1]+theta[2]))
h[1,2]=-n*(trigamma(theta[1]+theta[2]))
h[2,1]=-n*(trigamma(theta[1]+theta[2]))
h[2,2]=n*(trigamma(theta[2])-trigamma(theta[1]+theta[2]))
return(-h)
}
theta0<-c(3,60)
newton=nlminb(theta0,log.like,gradient,hessian,y=hr_average,n=length(hr_average))



Question 14.1
earth<-c(5.0,13.0,7.2,6.8,12.8,5.8,9.5,6.0,3.8,14.3,1.8,6.9,4.7,9.5)
clay<-c(0.9,12.9,2.6,3.5,26.6,1.5,13.0,8.8,19.5,2.5,9.0,13.1,3.6,6.9)
goodhue<-c(14.3,6.9,7.6,9.8,2.6,43.5,4.9,3.5,4.8,5.6,3.5,3.9,6.7)
basement1<-c(1,1,1,1,1,0,1,1,1,0,1,1,1,1)
basement2<-c(0,1,1,0,1,1,1,1,1,0,1,1,1,0)
basement3<-c(1,0,1,0,1,1,1,1,1,1,1,1,1)

