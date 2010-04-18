data <- read.table("../SAT.dat",header = TRUE)

logpost.mutau <- function(mu, tau, y=data$esteffect, sig=data$seeffect){
	logprior <- -log(tau)/2
	loglike  <- sum(dnorm(y, mean=mu, sd=sqrt(sig^2 + tau), log = TRUE))
	logpost <- logprior + loglike
	return(logpost)
}

vlogpost <- Vectorize(logpost.mutau)

n.grid <- 100
mu.grid <- ppoints(n.grid)*25-4.5
tau.grid <- ppoints(n.grid)*10


vals <- data.matrix(expand.grid(mu.grid, tau.grid))

logpost<-vlogpost(vals[,1],vals[,2])

post <- exp(logpost-max(logpost))
post <- post/sum(post)
summary(post)
mpost <- matrix(post, n.grid,n.grid)

contour(mu.grid,tau.grid, mpost, drawlabels=FALSE, main= expression(paste("Marginal Posterior Distribution of ", plain(p)(mu,tau^2))))