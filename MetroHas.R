require(gdata)
require(MASS)
library("matrixStats")

### we import the data
setwd("/Users/arnaudb/Desktop/montecarlo")
bos2 = read.xls ("Bostonhousing.xls", sheet=1, header=TRUE)
bos=bos2 = read.xls ("boston1.xls", sheet=1, header=TRUE)
bos <- subset(bos, select = -c(X,X.1,X.2,X.3,X.4) )
bos$cons <- matrix(1,nrow(bos),1)


### let's adjust data ###

Boston_copy<-cbind(data.matrix(bos2[,1:3]),data.matrix(bos2[,5:14]))

#Let's demean the covariates

for (j in 1:13) {
  m = mean(Boston_copy[,j])
  n = length(Boston_copy[,1])
  for (i in 1:n){
    Boston_copy[i,j] = Boston_copy[i,j] - m
  }
}

#Let's rescale the covariates
for (j in 1:13) {
  StandardDeviationSquared = 0
  n = length(Boston_copy[,1])
  for (i in 1:n){
    StandardDeviationSquared = StandardDeviationSquared + (Boston_copy[i,j]^2)/n
  }
  for (i in 1:n){	
    Boston_copy[i,j] = Boston_copy[i,j]/sqrt(StandardDeviationSquared)
  }
}

X<-Boston_copy[,1:12]
Y<-Boston_copy[,13]


#density multivariate normal distribution
densite <- function(x , beta , sig){
  d=length(beta)
  return((2*pi)^(-d) * (det(sig))^(-1/2) * exp(-( t(x-beta) * solve(sig) * (x-beta)) / 2) ) 
}


#compute the OLS, used as starting point
compute.beta <- function(X,Y){
  N=nrow(X)
  mat = matrix(0, ncol(X), ncol(X))
  vect= matrix (0, ncol(X), 1)
  for (k in 1:N){
    mat= mat + X[k,] %*% t(X[k,])
    vect= vect + Y[k] * as.matrix(X[k,])
  }
  
  return ( solve(mat) %*% vect)
}

compute.sig <- function(X,Y,sig2){
  N=nrow(X)
  mat = matrix(0, ncol(X), ncol(X))
  for (k in 1:N){
    mat= mat + unlist(X[k,]) %*% t(unlist(X[k,]))
  }
  return ( solve(1/sig2*mat) )
}

#a posteriori density (we drop all constants, they disappear when we calculate the fraction latter)
p <- function(beta, Y , X, sig2, lambda, kappa){
  N=length(Y)
  val= 0
  for (k in 1:N){
    val= val + (Y[k]- t(beta) %*% as.matrix(X[k,]))^2
  }
  #print(val/N)
  return ( exp(-lambda*sum(abs(beta)^kappa) - (val)/(2*sig2)))
}



metropolis_hastings <- function(X,Y,lambda,kappa,sig2,T =2000, scale){
  D=ncol(X)
  theta=matrix(0,D,T)
  theta[,1]=compute.beta(X,Y)
  SIGMA= scale*compute.sig(X,Y,sig2)
  accept=0
  for (k in 2:T){
    guess=mvrnorm(n = 1, theta[,k-1] , SIGMA)
    alpha=min(1, p(guess, Y, X, sig2 , lambda, kappa ) / p(theta[,k-1],  Y, X, sig2 , lambda, kappa ) )
    if (alpha == 1){
    }
    if (runif(1) < alpha){
      accept=accept+1
      theta[,k]=guess
    }else{
      theta[,k]=theta[,k-1]
    }
  }
  print(accept/T)
  return(theta)
}

theta <- metropolis_hastings(X,Y,1,1,1,20000,0.45) # acceptance rate: 0.246

## Let's choose the burn-in and check convergence
plot_cumsum <- function(theta,n){
  plot(cumsum(theta[n,])/seq_along(theta[n,]), type='l', xlab= "itérations", ylab="moyenne cumulée")
  grid()
}

par(mfrow=c(2,2))
for (k in 1:4){
  plot_cumsum(theta,k)
}
for (k in 5:8){
  plot_cumsum(theta,k)
}
for (k in 9:12){
  plot_cumsum(theta,k)
}


##We compute our final estimation of beta
beta_final <- function(theta, burnin){
  return(rowMeans(theta[,burnin:ncol(theta)])
)
}


beta_f=beta_final(theta, 15000)


#we compute the variance of the estimator
var_beta <- function(theta, burnin){
  return(rowVars(theta[,burnin:ncol(theta)]))
}


sqrt(var_beta(theta, 15000))




##Graph used as an example
ex <- metropolis_hastings(as.matrix(X[,12]),Y,1,1,1,2000,30)
plot_cumsum(ex, 1)
mean(ex[1200:2000])
compute.beta(as.matrix(X[,12]), Y)


##We plot the graph

plot(bos$LSTAT, bos$MV, pch = 10, cex = 0.5, col = "blue", xlab="", ylab="")
abline(34.55, -0.95, col='green') #OLS
abline(32.91, -0.858, col='red')  # LASSO
#For kappa=1: coef=-0.858 intercept: 32.91
#kappa=0: cof=-0.95 intercept=34.55



