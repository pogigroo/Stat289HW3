#######################################
R-Code for Problem-12
#######################################
Part-(b)
#######################################
meanX <- 30 
meanY <- 0  
varX  <- 5  
varY  <- 1  
rho   <- -0.4 # CORRELATION
covXY <- rho * sqrt(varX) * sqrt(varY)

mu <- c(meanX, meanY)
sigma <- matrix(c(varX, covXY, covXY, varY), 2, 2)
xgrid <- seq(20, 40, length = 100)
ygrid <- seq(-2, 2,  length = 100)
z <- matrix(NA, nrow = length(xgrid), ncol = length(ygrid))

for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- (1/(2*pi*det(sigma)))*exp( -1*((xgrid[i]-meanX)^2+(ygrid[j]-meanY)^2-2*rho*(xgrid[i]-meanX)*(ygrid[j]-meanY)) /(2*(1-rho^2)) )
  }
}
contour(x = xgrid, y = ygrid, z = z)

#######################################
Part-(e)
#######################################
y<-c(24,25,31,31,22,21,26,20,16,22)
t<-c(1,2,3,4,5,6,7,8,9,10)# same forcast will be found if we use complete years.
c<-lm(y~t)
summary(c)

#######################################
Part-(f)
#######################################

# function to estimate the Unnormalized density
fn <- function(a, b){
c1 <- exp(-a - b*ti)
c2<- (a + b*ti)^yi
c3<- factorial(yi)
out <- c1*c2/c3
prod(out)
}
# original data
ti <- c(1,2,3,4,5,6,7,8,9,10)
yi <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)
# contour plot on a grid of 1000*1000
l <- 1000   
a <- seq(18, 42, length = l)
b <-  seq(-3, .5, length = l)
# evaluation of fn on a grid of 1,000,000 values
z <- matrix(NA, ncol = l, nrow = l)
for(i in 1:l){
 for(j in 1:l){ z[i,j] <- fn(a[i], b[j])} 
      }
# the plot
contour(a, b, z, xlab = expression(alpha), ylab = expression(beta),
               xlim = c(18, 42), ylim = c(-3, .5), las = 1)
#######################################
Normalized density(nd).
finding the constant c and then multiplying it by z gives us normalized z.
we have to draw 1000 samples from marginal density of alpha
we have to draw 1000 samples from following nd.
#######################################
c<-1/sum(z)
nd<-c*z #normalized density for both alpha and beta.
dim(nd)
sum(nd)
alphaden<-colSums(nd) #marginal density for alpha. we have to draw alpha from this density
U<-runif(1000,min=0,max=1)
fa<-data.frame(a, alphaden, alphacdf = 0)

for (ndx in 2:l) {
  fa[ndx, "alphacdf"]<-fa[ndx-1, "alphacdf"]+fa[ndx, "alphaden"]
}
inversecdf<-function(u = .3, fa = fa) {
  return(fa$a[findInterval(u, fa[,3])])
}
alphaSample<-sapply(U, inversecdf, fa = fa)

    


#######################################
Part-(g)
#######################################





