data <- read.table("../SAT.dat",header = TRUE)

logpost.mutau <- function(mu, tau, y=data$esteffect, sig=data$seeffect){
	logprior <- -log(tau)/2
	loglike  <- sum(dnorm(y, mean=mu, sd=sqrt(sig^2 + tau), log = TRUE))
	logpost <- logprior + loglike
	return(logpost)
}

mugrid <- ppoints(1000)*20-3
taugrid <- ppoints(1000)*10 #per the hint

logpost <- matrix(NA,1000,1000)
for (i in 1:1000){ 
  for (j in 1:1000){
    logpost[i,j] <- logpost.mutau(mugrid[i],taugrid[j])
  }
}

post <- exp(logpost-max(logpost))
post <- post/sum(post)

par(mar=c(3.5,3.5,2,1), mgp=c(2,.65,0), las=1, cex.main=1.5)
contour(mugrid,taugrid,post, ylim=c(0,0.6), drawlabels=FALSE, main= expression(paste("Marginal Posterior Distribution of ", plain(p)(mu,tau^2))))