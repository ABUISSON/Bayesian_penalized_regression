#install.packages("BayesBridge")
library("BayesBridge")


Boston <- read.csv2("C:/Users/bessem/Documents/ENSAE/2A/S2/Monte Carlo/Boston housing2.csv")

library("Matrix")

Boston_copy<-cbind(data.matrix(Boston[,1:3]),data.matrix(Boston[,5:14]))

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


#Initialisation
lambda = 1
sigma2=1

X<-Boston_copy[,1:12]
Y<-Boston_copy[,13]

modele<- bridge.reg(Y, X, 50000, alpha=1,
           sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,
           alpha.a=1.0, alpha.b=1.0,
           sig2.true=1.0, tau.true=10,
           burn=1000, method="stable", ortho=FALSE)

summary(modele$beta)

beta <-modele$beta
beta[10000,1:12]
plot(9900:10000,beta[9900:10000,1], "l",grid())

boxplot(beta[,1],outline=FALSE, border="blue")