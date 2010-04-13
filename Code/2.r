data <- read.table("../SAT.dat",header = TRUE)

logpost.mutau <- function(mu, tau, y=data$esteffect, sig=data$seeffect){
	logprior <- -log(tau)/2
	loglike  <- sum(log(dnorm(y, mean=mu, sd=sqrt(sig^2 + tau))))
	logpost <- logprior + loglike
	return(logpost)
}

mugrid <- ppoints(100)*20-5
taugrid <- ppoints(100)*10 #per the hint

logpost <- matrix(NA,100,100)
for (i in 1:100){ 
  for (j in 1:100){
    logpost[i,j] <- logpost.mutau(mugrid[i],taugrid[j])
  }
}

post <- exp(logpost-max(logpost))
post <- post/sum(post)

contour(mugrid,taugrid,post, drawlabels=FALSE)