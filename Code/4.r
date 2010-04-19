data <- read.table("../SAT.dat",header = TRUE)
logpost.mutau <- function(mu, tau, y=data$esteffect, sig=data$seeffect){
	logprior <- -log(tau)/2
	loglike  <- sum(dnorm(y, mean=mu, sd=sqrt(sig^2 + tau), log = TRUE))
	logpost <- logprior + loglike
	return(logpost)
}

vlogpost <- Vectorize(logpost.mutau)

n.grid <- 500
mu.grid <- ppoints(n.grid)*25-4.5
tau.grid <- ppoints(n.grid)*10

vals <- data.matrix(expand.grid(mu.grid, tau.grid))

logpost<-vlogpost(vals[,1],vals[,2])
post <- exp(logpost-max(logpost))
post <- post/sum(post)

# ind <- 1:length(post)
samp.ind <- sample.int(length(post), size = 1000, replace = TRUE, prob = post)
samp.post <- vals[samp.ind, ]

rpost.theta <- function(mu=samp.post[,1], tausq=samp.post[,2], y, sig, size = length(mu)){
	V <- 1 / (1/tausq+1/sig^2)
	thetahat <- (mu/tausq+y/sig^2)*V
	draws <- rnorm(size, thetahat, sqrt(V))
	return(draws)
}

theta.draws <- mapply(rpost.theta, y=data$esteffect, sig=data$seeffect, SIMPLIFY = FALSE)

mean(data$esteffect)
rbind(sapply(theta.draws, mean), y=data$esteffect, sigsq=data$seeffect)
