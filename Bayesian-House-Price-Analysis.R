library(mlbench)
library(R2OpenBUGS)
library(coda)

# Housing prices are a trendy news topic right now, and Bayesian statistics is a trendy field 
# in stats, so let's put them together and be double trendy.

WorkDir<- "C:\\Users\\Viktor\\Desktop"
setwd(WorkDir)
set.seed(42)

# load our dataset
data(BostonHousing)
dat = BostonHousing

J <- nrow(dat)
y <- dat$medv # median value of owner-occupied homes in USD 1000's
crim <- dat$crim # per capita crime rate by town
zn <- dat$zn # 	proportion of residential land zoned for lots over 25,000 sq.ft
indus <- dat$indus # proportion of non-retail business acres per town
chas <- as.numeric(dat$chas) # Charles River dummy variable
nox <- dat$nox # nitric oxides concentration (parts per 10 million)
rm <- dat$rm # average number of rooms per dwelling
age <- dat$age # proportion of owner-occupied units built prior to 1940
dis <- dat$dis # weighted distances to five Boston employment centres
rad <- dat$rad # index of accessibility to radial highways
tax <- dat$tax # full-value property-tax rate per USD 10,000
ptratio <- dat$ptratio # pupil-teacher ratio by town
b <- dat$b # 1000(B - 0.63)^2 where B is the proportion of town population who identify as black
lstat <- dat$lstat # 	percentage of lower status of the population

bug.dat <- list("J", "y", "crim", "zn", "indus", "chas", "nox", "rm", "age", "dis", "rad",    
                "tax", "ptratio", "b", "lstat")


cat("
    model {
    for (j in 1:J)
    {
    y[j] ~ dnorm(mu[j], tau)
    mu[j] <- B[1]*crim[j] + B[2]*zn[j] + B[3]*indus[j] + B[4]*chas[j] +
    B[5]*nox[j] + B[6]*rm[j] + B[7]*age[j] + B[8]*dis[j] + B[9]*rad[j] + 
    B[10]*tax[j] + B[11]*ptratio[j] + B[12]*b[j] + B[13]*lstat[j] + B[14]
    }
    
    tau~dgamma(.5,.01)
    for(i in 1:14){
    B[i]~dnorm(0,0.0625)}
    }", file="BostonModel.txt")


inits<-function(){ list(B=rnorm(14), tau=runif(.5,1))} 
params=c("B", "tau")


housing.sim <- bugs(bug.dat, inits, model.file = "C:\\Users\\Viktor\\Desktop\\BostonModel.txt",
                    params, n.chains=3, n.iter=20000, n.burnin=5000, n.thin=3)

print(housing.sim, digits.summary = 3)

SArray = housing.sim$sims.array


par(mfrow=c(1,1))
#acf plots show correlation quickly converges to 0
acf(SArray[,1,"B[1]"], main="Sample ACF Plot")

#trace plots show convergance as well, hovering nicely around a central line
matplot(1:15000,SArray[,,"B[1]"], main="Sample Trace Plot",xlab="index",type="l")

geweke.diag(SArray[,,"B[1]"], frac1=0.1, frac2=0.5)
# with large enough number of iterations, T test approximated by standard normal
# absolute value of z scores for all chains in all variables is less than 1.96
# this means we can treat the first 10% of the data as burn-in

housingcoda.sim <- bugs(bug.dat, inits, model.file = "C:\\Users\\Viktor\\Desktop\\BostonModel.txt",
                        params, n.chains=3, n.iter=20000, n.burnin=5000, n.thin=3, codaPkg = TRUE)
housing.coda <- read.bugs(housingcoda.sim)
gelman.diag(housing.coda, confidence = 0.95, transform=FALSE, autoburnin=FALSE,
            multivariate=TRUE) # shrink factors all very close to 1.

gelman.plot(housing.coda) # visibly also seems to have converged



# posterior distribution of our Beta parameters

par(mfrow=c(2,2))

plot(c(-0.5,0.5),c(0,15), type="n",main="crim")
apply(SArray[,,"B[1]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

plot(c(-0.5,0.5),c(0,35), type="n",main="zn")
apply(SArray[,,"B[2]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

plot(c(-0.5,0.5),c(0,15), type="n",main="indus")
apply(SArray[,,"B[3]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

plot(c(-0.5,8),c(0,1), type="n",main="chas")
apply(SArray[,,"B[4]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")


####

plot(c(-16,6),c(0,0.5), type="n",main="nox")
apply(SArray[,,"B[5]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

plot(c(0,7),c(0,2), type="n",main="rm")
apply(SArray[,,"B[6]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

plot(c(-0.5,0.5),c(0,35), type="n",main="age")
apply(SArray[,,"B[7]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

plot(c(-2.5,0),c(0,3), type="n",main="dis")
apply(SArray[,,"B[8]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

####

plot(c(-0.2,0.7),c(0,15), type="n",main="rad")
apply(SArray[,,"B[9]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

plot(c(-0.2,0.2),c(0,120), type="n",main="tax")
apply(SArray[,,"B[10]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

plot(c(-1.5,0),c(0,6), type="n",main="ptratio")
apply(SArray[,,"B[11]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

plot(c(-0.1,0.1),c(0,160), type="n",main="b")
apply(SArray[,,"B[12]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")

###

plot(c(-1,0.5),c(0,10), type="n",main="lstat")
apply(SArray[,,"B[13]"], 2, function(x)lines(density(x)))
abline(v = 0, col="red")
