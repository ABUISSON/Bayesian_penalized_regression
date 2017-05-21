######### Prépartion des données ########

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

X<-Boston_copy[,1:12]
Y<-Boston_copy[,13]

##### Densité a priori, vraisemblance et densité a posteriori #####

a_priori <- function(beta, lambda, kappa){
	p = length(beta)
	S = 0
	for (i in 1:p){
		S = S + (abs(beta[i])^kappa)
		}
	x = exp(-lambda*S)
	return(x)
	}

vraisemblance <- function(Y,X,beta,sigma){
	N= length(X[,1])
	const = (N/2)*log(1/(2*pi*sigma))
	S = 0
	for (i in 1:N){
		S = S + ((Y[i]- t(beta)%*%X[i,])^2)
		}
	v = exp(const - (1/(2*(sigma^2)))*S)
	return(v)
	}

a_posteriori <- function(beta,Y,X,lambda,kappa,sigma){
	A = a_priori(beta, lambda, kappa)
	B = vraisemblance(Y,X,beta,sigma)
	C = A*B
	return(C)
	}

##### Simulons M beta avec kappa = 0 #####

library("mvtnorm")

Simulation <- function(M,Y,X,lambda,sigma){
	N = length(X[,1])
	S = 0
	S2 = 0
	for (i in 1:N){
	S = S + Y[i]*X[i,]
	S2 = S2 + X[i,]%*%t(X[i,])  
		}
	S2 = S2 * (1/(sigma^2))
	Var = solve(S2)
	beta_average = Var %*% S
	beta_simu = rmvnorm(M, beta_average, Var)
	return(beta_simu/(sigma^2))
	} 

#On simule un grand nombre de beta avec kappa = 0
lambda = 1
sigma = 1
kappa=1
beta_simu = Simulation(50000,Y,X,lambda,sigma)


###hist(beta_simu[,1]) ##### avec kappa = 0

######## Importance sampling (après simplification des calculs) ########

Importance_Sampling_V2 <- function(beta_simulation,lambda,kappa){
	MC = 0
	Norm = 0
	M1 = length(beta_simulation[,1]) 
	for (i in 1:M1){
		MC = MC + beta_simu[i,]*a_priori(beta_simulation[i,], lambda, kappa)
		Norm = Norm + a_priori(beta_simulation[i,], lambda, kappa)
	}
	beta_average = MC / Norm 
	return(beta_average)
	}

Erreur_Importance_Sampling_V2 <- function(beta_simulation,lambda,kappa){
	beta_average = Importance_Sampling_V2(beta_simulation,lambda,kappa)
	N = 0
	Norm = 0 
	M = length(beta_simulation[,1])
	for (i in 1:M){
		N = N + (a_priori(beta_simulation[i,], lambda, kappa)*(beta_simulation[i,]-beta_average))**2
		Norm = Norm + a_priori(beta_simulation[i,], lambda, kappa)
		}
	Norm_squared = Norm**2
	N = N*M/Norm_squared
	return(N**0.5)
	}

##### Intervalle de confiance asymptotique #####
Intervalle_de_confiance_95 <- function(beta_simulation,lambda,kappa){
	beta_average = Importance_Sampling_V2(beta_simulation,lambda,kappa)
	Erreur = Erreur_Importance_Sampling_V2(beta_simulation,lambda,kappa)
	M = length(beta_simulation[,1])
	Intervalle <- cbind(beta_average,beta_average) ## crée un tableau de taille 12*2
	for (i in 1:12){
		Intervalle[i,1] = beta_average[i]-1.96*Erreur[i]/(M**0.5)
		Intervalle[i,2] = beta_average[i]+1.96*Erreur[i]/(M**0.5)
		}
	return(Intervalle)
	}

####### tests importance sampling ####### 

lambda = 1
sigma = 1
kappa=1

b2 = Importance_Sampling_V2(beta_simu,lambda,kappa)
c2 = Erreur_Importance_Sampling_V2(beta_simu,lambda,kappa)
Intervalle = Intervalle_de_confiance_95(beta_simu,lambda,kappa)

print(b2)
print(c2)
##print(Intervalle) ##intervalle de confiance à 95%
Tableau_importance <- cbind(Intervalle[,1],b2,Intervalle[,2])
print(Tableau_importance)

######## Effective Sample Size ########

ESS <- function(beta_simulation,lambda,kappa){
  Numerator = 0
  Denominator = 0
  M = length(beta_simulation[,1])
  for (i in 1:M){
    Numerator = Numerator + a_priori(beta_simulation[i,], lambda, kappa)
    Denominator = Denominator + a_priori(beta_simulation[i,], lambda, kappa)**2            
  }
  Numerator_squared = Numerator**2
  return(Numerator_squared/Denominator)
}

##### test pour les lambda #####

###beta_simu = Simulation(50000,Y,X,1,sigma)
ESS1 = ESS(beta_simu,1,kappa)
print(ESS1)

######## Resampling ########

###### Calculs du poids à attribuer à chaque beta simulé ######

weights <- function(beta_simulation,lambda,kappa){
	weight_total = a_priori(beta_simulation[1,], lambda, kappa)
	weights_beta <- c(weight_total)
	M = length(beta_simulation[,1])
	for (i in 2:M){
		weight_new = a_priori(beta_simulation[i,], lambda, kappa)
		weight_total = weight_total + weight_new
		weights_beta <- cbind(weights_beta ,weight_new)
		}
	for (i in 1:M){
		weights_beta[i] =  weights_beta[i]/weight_total
		}
	return(weights_beta)
	} 

##W = weights(beta_simu,lambda,kappa)
##N = length(W)
##print(length(W))

#install.packages("multinomRob")
library("multinomRob")

##s = rmultinomial(N,W)
#print(s)

##### créons les nouveaux beta suivants la loi a posteriori avec kappa #####

Sampling <- function(beta_simulation,lambda,kappa){ 
	W = weights(beta_simulation,lambda,kappa)
	N = length(W)
	s = rmultinomial(N,W)
	a = 1
	beta_sampled <- beta_simulation
	for (i in 1:N){
		if (s[i] != 0){
			for (j in a:(a+s[i]-1)){
				beta_sampled[j,]<-beta_simulation[i,]
				}
			a = a + s[i]
			}
		}
	return(beta_sampled)
	}

##beta_sampled = Sampling(beta_simu,lambda,kappa) 

##### Analyse de la moyenne et de la variance après resampling ######

Average_resampling <- function(beta_kappa){
	MC = 0
	M1 = length(beta_kappa[,1]) 
	for (i in 1:M1){
		MC = MC + beta_kappa[i,]
	}
	beta_average = MC / M1
	return(beta_average)
	}

Erreur_resampling <- function(beta_kappa){
	beta_average = Average_resampling(beta_kappa)
	N = 0 
	M = length(beta_kappa[,1])
	for (i in 1:M){
		N = N + (beta_kappa[i,]-beta_average)**2
		}
	N = N/M
	return(N**0.5)
	}

##### Intervalle de confiance asymptotique #####
Intervalle_de_confiance_95_resampling <- function(beta_kappa,lambda,kappa){
	beta_average = Average_resampling(beta_kappa)
	Erreur = Erreur_resampling(beta_kappa)
	M = length(beta_kappa[,1])
	Intervalle <- cbind(beta_average,beta_average) ## crée un tableau de taille 12*2
	for (i in 1:12){
		Intervalle[i,1] = beta_average[i]-1.96*Erreur[i]/(M**0.5)
		Intervalle[i,2] = beta_average[i]+1.96*Erreur[i]/(M**0.5)
		}
	return(Intervalle)
	}


###### test resampling ######
lambda = 1
##beta_simu = Simulation(50000,Y,X,lambda,sigma)

beta_sampled = Sampling(beta_simu,lambda,kappa)
b3 = Average_resampling(beta_sampled)
c3 = Erreur_resampling(beta_sampled)

print(b3)
print(c3)
Intervalle_resampling = Intervalle_de_confiance_95_resampling(beta_sampled,lambda,kappa)
##print(Intervalle_resampling)
Tableau_resampling <- cbind(Intervalle_resampling[,1],b3,Intervalle_resampling[,2])
print(Tableau_resampling)


#### exploitation du Resampling ####
##Densités a posteriori univariées des coefficients des variables 1,6,11 et 12 pour Kappa =1

par(mfrow=c(2,2))
hist(beta_sampled[,1], main = "", 
	ylab="effectif", xlab="densité",breaks = 50)

hist(beta_sampled[,6], main = "", 
     ylab="effectif", xlab="densité",breaks = 50)

hist(beta_sampled[,11], main = "", 
     ylab="effectif", xlab="densité",breaks = 50)


hist(beta_sampled[,12], main = "", 
     ylab="effectif", xlab="densité",breaks = 50)

