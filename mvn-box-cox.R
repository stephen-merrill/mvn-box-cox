dat <- read.table('https://tofu.byu.edu/stat666/datasets/oliver3b.txt', header = TRUE)
library(MASS)

#EDA
  dim(dat) #206 observations
  apply(dat,2, mean)
  apply(dat,2, median)
  #boxplots
  boxplot(dat[,c(1,5)])
  boxplot(dat[,c(2,3,6,7,8)])
  boxplot(dat[,4])

################################################
#How to check for Normality
################################################
#lambda 
lam <- function(x, lambda){
  lambs <- matrix(0, ncol=length(lambda), nrow=length(x))
  for(i in 1:length(lambda)){
    if(lambda[i] != 0){
      lambs[,i] <- (x^(lambda[i])-1)/(lambda[i])
    }else{
      lambs[,i] <-log(x)
    }
  }
  return(lambs)
}
lam(dat[,1], c(0,1,2,3,4))

#Box-Cox likelihood
boxy <- function(x, lambda){
  lambdas <- lam(x, lambda)
  like <- numeric(length(lambda))
  for(i in 1:ncol(lambdas)){
    s2lam <- (1/length(x))*sum((lambdas[,i]-mean(lambdas[,i]))^2)
    like[i] <- -length(x)/2*log(s2lam) + (lambda[i]-1)*sum(log(x))
  }
  return(like)
}
seqs <- seq(-10, 10, by = 0.01)
evals <- boxy(dat[,1], seqs)
plot(seqs, evals, type = 'l')
seqs[which.max(evals)]

evals <- numeric(8)
for(i in 1:8){
  temp <- boxy(dat[,i], seqs)
  evals[i] <- seqs[which.max(temp)]
}
first <- evals

############################################
#MULTIVARIATE BOXCOX
############################################
#evals are lambdas
multivbc <- function(evals) {
  n <- 206
  evals2 <- numeric(8)
  MLE <- matrix(0, ncol = 3, nrow = 8)
  
  for(i in 1:length(evals)){
    eigs <- zapsmall(c(evals[i], evals[i] + .001, evals[i] -0.001))
    dats <- matrix(0, ncol =8, nrow = 206)
    for(j in 1:8){
      dats[,j]<- lam(dat[,j], evals[j])
    }
    #test which of the 3 lambda is best
    for(k in 1:3){
      evals[i] <- eigs[k]
      dats[,i] <- lam(dat[,i],eigs[k])
      #estimate S with cov
      #Plug into MLE
     MLE[i,k] <- -n/2 * log(det(cov(dats))) + sum((evals-1)*apply(dat,2, function(x) sum(log(x))))
    }
    evals2[i] <- eigs[which.max(MLE[i,])]
  }
  evals2
}
##################################################
#Iterate until all 8 rows don't change
comp <- rep(FALSE,8)
iter <- 0
while(!all(comp)) {
  print(first)
  current <- multivbc(first)
  print(current)
  first <- multivbc(current)
  comp <- first==current
  iter <- iter+1
}
###############################################
#Normality Check
###############################################
par(mfrow=c(2,4))
for(i in 1:8) {
  hist(dat[,i]^evals[i],20,main=colnames(dat)[i])
}
for(i in 1:8) {
  hist(dat[,i]^first[i],20,main=colnames(dat)[i])
}

###############################################
#How to identify outliers
###############################################
dat.trans <- dat
for(i in 1:ncol(dat)) {
  dat.trans[,i] <- (dat[,i]^first[i]-1)/first[i]
}
D2 <- numeric(nrow(dat))
Sinv <- solve(cov(dat.trans))
for(i in 1:nrow(dat)) {
  D2[i] <- t(t(dat.trans[i,]-colMeans(dat.trans)))%*%Sinv%*%t(dat.trans[i,]-colMeans(dat.trans))
}
D2[which(D2>20)]
n <- nrow(dat)
p <- ncol(dat)
a <- p/2 
b <- (n-p-1)/2
alpha <- (p-2)/(2*p)
beta <- (n-p-3)/(2*(n-p-1))
v <- qbeta(((1:n) - alpha)/(n-alpha-beta+1),a,b)
u <- (n*D2)/((n-1)^2)
which(u>0.12)
par(mfrow=c(1,1))
plot(v,sort(u),ylab=expression(paste("Scaled ", D^2)),xlab="Quantiles")
abline(0,1)
plot(v[-206],sort(u)[-206])
plot(v[-c(204:206)],sort(u)[-c(204:206)])
