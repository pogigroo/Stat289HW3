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
# hist(samp.ind)
samp.post <- vals[samp.ind, ]

summary(samp.post[,1])
mean(samp.post[,1])
quantile(samp.post[,1],c(0.025,0.5,0.975))

summary(samp.post[,2])
mean(samp.post[,2])
quantile(samp.post[,2],c(0.025,0.5,0.975))

#points(quantile(vals[,1],c(0.025,0.5,0.975)), quantile(vals[,2],c(0.025,0.5,0.975)), pch=19)

# boxplot(list(vals[,1], vals[,2]), names=c("mu","tausq"))

