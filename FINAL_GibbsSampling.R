# Bayesian Lasso Gibbs sampling
library("Matrix")
library("monomvn")

#donnees Boston housing dataset
Boston <- read.csv2("/Users/rimanetoumi/Desktop/Projet MC/boston housing2.csv", sep = ";")
Boston_copy<-cbind(data.matrix(Boston[,1:3]),data.matrix(Boston[,5:14]))

#Centrer les variables explicatives
for (j in 1:13) {
  m = mean(Boston_copy[,j])
  n = length(Boston_copy[,1])
  for (i in 1:n){
    Boston_copy[i,j] = Boston_copy[i,j] - m
  }
}

#Normalisation de la variance des variables explicatives
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

n<-1000

#Gibbs Sampling
modele <- blasso(X,Y,T=n,RJ=FALSE,beta=NULL,lambda2 = 1, s2 = 1, mprior=0, rd = FALSE, ab=NULL,theta=0, rao.s2=FALSE, icept=FALSE, normalize=FALSE, verb=1)
summary(modele)
beta <-modele$beta

#Plot des coefficients
plot(900:1000,beta[900:1000,1], "l")
plot(900:1000,beta[900:1000,2], "l")
plot(900:1000,beta[900:1000,3], "l")
plot(900:1000,beta[900:1000,4], "l")
plot(900:1000,beta[900:1000,5], "l")
plot(900:1000,beta[900:1000,6], "l")
plot(900:1000,beta[900:1000,7], "l")
plot(900:1000,beta[900:1000,8], "l")
plot(900:1000,beta[900:1000,9], "l")
plot(900:1000,beta[900:1000,10], "l")
plot(900:1000,beta[900:1000,11], "l")
plot(900:1000,beta[900:1000,12], "l")


#Convergence des estimateurs 
    #plot de la moyenne cumulée 
plot_cumsum <- function(theta,n){
  plot(cumsum(theta[n,])/seq_along(theta[n,]), type='l', xlab= "itérations", ylab="moyenne cumulée")
  grid()
}

plot_cumsum(t(beta),1)
plot_cumsum(t(beta),2)
plot_cumsum(t(beta),3)
plot_cumsum(t(beta),4)
par(mfrow=c(2,2)) #pour plotter 4 par 4
plot_cumsum(t(beta),5)
plot_cumsum(t(beta),6)
plot_cumsum(t(beta),7)
plot_cumsum(t(beta),8)
par(mfrow=c(2,2)) 
plot_cumsum(t(beta),9)
plot_cumsum(t(beta),10)
plot_cumsum(t(beta),11)
plot_cumsum(t(beta),12)
par(mfrow=c(2,2)) 


#Moyenne, écart type, intervalles de confiance à 95% des beta estimés
for (i in 1:12) {
  print(i)
  print("moyenne : ")
  print(mean(beta[400:1000 ,i]))
  print("écart type :")
  print(sqrt(var(beta[400:1000,i])))
  print("intervalle de confiance à 95% : ")
  print(quantile(beta[400:1000,i], c(0.05, 0.95)))
}




