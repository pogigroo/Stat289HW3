# Simulate x~exp(lambda) using inverse-CDF method
lambda = 3
set.seed = 100
u = runif(1000)
x = -log(u)/lambda

# Theoretical and empirical histogram
postscript("lecture1-exp.ps",horiz=T,height=6,width=6)
par(mar=c(3,3,1,1),las=1)
hist(x,prob=T)
xg = seq(0,4,len=200)
lines(xg,dexp(xg,lambda),col=2)
dev.off()

# Theoretical and empirical quantiles
quant.theor = round(qexp(c(.025,.50,.975),lambda),3)
quant.empir = round(quantile(x,c(.025,.50,.975)),3)

# Theoretical and empirical moments
moment.theor = rep(round(1/lambda,3),2)
moment.empir = round(c(mean(x),sd(x)),3)

cat("\n")
cat("Theoretical Quantiles:", quant.theor,"\n")
cat("Empirical Quantiles:  ", quant.empir,"\n\n")
cat("Theoretical Moments:", moment.theor,"\n")
cat("Empirical Moments:  ", moment.empir,"\n\n")
