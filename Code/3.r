data <- read.table("../SAT.dat",header = TRUE)

logpost.mutau <- function(mu, tau, y=data$esteffect, sig=data$seeffect){
	logprior <- -log(tau)/2
	loglike  <- sum(dnorm(y, mean=mu, sd=sqrt(sig^2 + tau), log = TRUE))
	logpost <- logprior + loglike
	return(logpost)
}

mugrid <- ppoints(1000)*20-3
taugrid <- ppoints(1000)*4 #smaller, to aid numerical integration

logpost <- matrix(NA,1000,1000)
for (i in 1:1000){ 
  for (j in 1:1000){
    logpost[i,j] <- logpost.mutau(mugrid[i],taugrid[j])
  }
}

post <- exp(logpost-max(logpost))
post <- post/sum(post)

######################################################################
# Numerical integration is easier than deriving the marginal 
# conditionals
######################################################################
post.mu.weights <- apply(post,1,sum)

post.mu.draws <- sample(mugrid, size=1000, replace=T, prob=post.mu.weights)
post.tau.draws <- array(NA, 1000)
for (i in post.mu.draws) {
	post.tau.draws[i] <- sample(taugrid, size=1, prob=post[mugrid==i, ])
}

mean(post.mu.draws)
quantile(post.mu.draws,c(0.025,0.5,0.975))
mean(post.tau.draws)
quantile(post.tau.draws,c(0.025,0.5,0.975))
summary(post.tau.draws)
hist(post.tau.draws)

######################################################################
# Numerically integrating over mu isn't as stable as the above
######################################################################
# post.tau.weights <- apply(post,2,sum)
# post.tau.draws <- sample(taugrid, size=1000, replace=T, prob=post.tau.weights)
# post.mu.draws <- array(NA, 1000)
# for (i in post.tau.draws) {
# 	post.mu.draws[i] <- sample(mugrid, size=1, prob=post[,taugrid==i])
# }
# 
# summary(post.mu.draws)
# summary(post.tau.draws)