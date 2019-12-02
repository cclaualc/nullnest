library('RSpectra')
library(MASS)


#############PROGRAM Spectral Radius###########
#This program calculates the following three quantities:
      #i) The spectral radius of a real network
      #ii) Its average null expectation using the theoretical expression derived in [1]
      #iii) The standard deviation of the distribution of spectral radius in the null ensemble, calculated analytically as derived in [1]
#The program takes as input the ordered real bi-adjacency matrix as well as the probability matrix over the null ensemble
#Code made by Clàudia Payrató Borràs. For any questions on the code, you may contact me by mail: claudiapb13@gmail.com
#When using this program, please cite: 
#[1] "Breaking the Spell of Nestedness: "The Entropic Origin of Nestedness in Mutualistic Systems", C. Payrató-Borràs et al. Phys. Rev. X 9 (2019)


########################################################
## Calculation of the REAL spectral radius #############
########################################################

##Reading input argument: index of the network to study
args <- commandArgs(TRUE)
intw <- as.integer(args[1])

#File with the real network
file_real = paste("matrix",c(intw),"ord.txt",sep="")

#Construct the real adjacency matrix
m = read.table(file_real)
m = as.matrix(m)
M = mat.or.vec(nrow(m)+ncol(m),nrow(m)+ncol(m))
M[1:nrow(m),(nrow(m)+1):(nrow(m)+ncol(m))] = m
M[(nrow(m)+1):(nrow(m)+ncol(m)),1:nrow(m)] = t(m)

#Calculate the real spectral radius
eig = abs(eigs_sym(M,1,which = "LM",opt=list(tol=1e-50,maxitr=100000))$values)



##########################################################
## Analytic calculation of the AVERAGE spectral radius ##
##########################################################

#File with the probability matrix of interactions
file_rand = paste("matrix",c(intw),"rand.txt",sep="")

#Construction of the adjacency matrix
m = read.table(file_rand)
m = as.matrix(m)
M = mat.or.vec(nrow(m)+ncol(m),nrow(m)+ncol(m))
M[1:nrow(m),(nrow(m)+1):(nrow(m)+ncol(m))] = m
M[(nrow(m)+1):(nrow(m)+ncol(m)),1:nrow(m)] = t(m)

#Calculation of the average spectral radius
eigav = abs(eigs_sym(M,1,which = "LM",opt=list(tol=1e-50,maxitr=100000))$values)



######################################
## Standard deviation calculation ##
######################################
#Here we calculate the standard deviation of the spectral radius using the expression obtained in [1]

#Construction of the Q-matrix
Q = eigav*diag(nrow(M)) - M

#Calculation of the Moore-Penrose inverse of Q
Qinv = ginv(Q,tol = sqrt(.Machine$double.eps))

#Calculation of the derivative matrix D
Qmult=Q%*%Qinv
D=t(diag(nrow(M))-Qmult)

#Calculation of the standard deviation using the expression obtained in [1]
i = 1
sigma = 0
while (i <= nrow(D)){
j = 1
while (j <= ncol(D)){
d = D[i,j]
sig = M[i,j]*(1.0-M[i,j])
sigma = sigma + d*d*sig
j = j + 1
}
i = i + 1
}
stdev = sqrt(sigma)



################
#Writing results
################

results = "analytic_spectral_radius.txt"
real_spectral_radius = eig
average_spectral_radius = eigav
standard_deviation = stdev
R = mat.or.vec(1,3)
R = rbind(real_spectral_radius, average_spectral_radius, standard_deviation)
write.matrix(t(R), file=results, sep=" ")



