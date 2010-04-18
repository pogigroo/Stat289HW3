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

post.mu.weights <- apply(post,1,sum)
post.tau.weights <- apply(post,2,sum)

# plot(post.mu.weights, type = "l")
# plot(post.tau.weights, type="l")
# 
# post.mu.draws <- sample(mugrid, size=1000, replace=T, prob=post.mu.weights)
# post.tau.draws <- array(NA, 1000)
# for (i in post.mu.draws) {
# 	post.tau.draws[i] <- sample(taugrid, size=1, prob=post[mugrid==i, ])
# }

post.tau.draws <- sample(taugrid, size=1000, replace=T, prob=post.tau.weights)
post.mu.draws <- array(NA, 1000)
for (i in post.tau.draws) {
	post.mu.draws[i] <- sample(mugrid, size=1, prob=post[,taugrid==i])
}

summary(post.mu.draws)
summary(post.tau.draws)

# par(mar=c(3.5,3.5,2,1), mgp=c(2,.65,0), las=1, cex.main=1.5)
# contour(mugrid,taugrid,post, ylim=c(0,0.6), drawlabels=FALSE, main= expression(paste("Marginal Posterior Distribution of ", plain(p)(mu,tau^2))))