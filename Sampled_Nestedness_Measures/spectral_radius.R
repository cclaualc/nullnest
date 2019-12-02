library('rARPACK')
library(MASS)
library('bipartite')

##########################PROGRAM Spectral Radius########################
#This program calculates the nestedness metrics Spectral Radius
#It produces the following information:
      #i) The spectral radius of a real network provided by the user
      #ii) The average spectral radius over a sampling of a null ensemble
      #iii) The standard deviation of the distribution of the spectral radius over the ensemble
##The program takes as input the ordered real bi-adjacency matrix as well as the probability matrix over the null ensemble
#Code made by Clàudia Payrató Borràs. For any questions on the code, you may contact me by mail: claudiapb13@gmail.com
#When using this program, please cite: 
#[1] "Lost in Nestedness? Assessing the performance of different metrics" Payrató-Borràs et al., 2019.


##Input arguments
args <- commandArgs(TRUE)
#Index 'intw', which identifies the network to study
intw <- as.integer(args[1])
#Size of the sampling: number of null networks to generate. We recommend a minimum size of 10³ networks
N_sampl <- as.integer(args[2])


########################################################
########Calculation of the REAL Spectral Radius##########
########################################################

#Opening the file containing the probability matrix
file_real = paste("matrix",intw,"ord.txt",sep="")

#Reading the real bi-adjacency matrix m
m = read.table(file_real)
m = as.matrix(m)

#Constructing the real adjacency matrix M
M = mat.or.vec(nrow(m)+ncol(m),nrow(m)+ncol(m))
M[1:nrow(m),(nrow(m)+1):(nrow(m)+ncol(m))] = m
M[(nrow(m)+1):(nrow(m)+ncol(m)),1:nrow(m)] = t(m)

#Calculating the spectral radius of M
real_spectral_radius = abs(eigs_sym(M,1,opt=list(tol=1e-50,maxitr=10000))$values)


###############################################################################
########Calculation of the AVERAGE Spectral radius over the ensemble##########
##############################################################################

#Opening the file containing the probability matrix
file_rand = paste("matrix",intw,"rand.txt",sep="")

#Reading the probability matrix of interactions, p
p = read.table(file_rand)

#Initializing vectors and parameters
set.seed(123456) #initializing the random number generator
toteig = 0.0 #sum over the spectral radius of each null network
toteigsq = 0.0 #sum over the squared spectral radius of each null network

#Start of the iterating loop to produce the sample
iter = 0 #index for the sampling
while (iter < N_sampl){

#Generation of a new null network b, using the matrix probability p
b = runif(nrow(p)*ncol(p))
b = matrix(b,nrow(p),ncol(p))
aux = b>p
b[aux]=0
b[!aux] = 1
b = as.matrix(b)
b = empty(b, count=FALSE) #Removing nodes with zero-degree (no interactions)s

#Construction of the adjacency matrix B
B = mat.or.vec(nrow(b)+ncol(b),nrow(b)+ncol(b))
B[1:nrow(b),(nrow(b)+1):(nrow(b)+ncol(b))] = b
B[(nrow(b)+1):(nrow(b)+ncol(b)),1:nrow(b)] = t(b)

#Calculation of the Spectral Radius of the null network
log = 1
eigone = tryCatch({
#try part
eigs_sym(B,1,which="LM",opt=list(tol=1e-50,maxitr=10000))$values },
#error part
error = function(e){
log = 0
},
warning = function(w){
log = 0
})

if (log == 1){ 
toteig = toteig + abs(eigone/N_sampl)
toteigsq = toteigsq + abs(eigone*eigone/N_sampl)
iter = iter + 1
}

} #end of the sampling loop

#Final calculations
average_spectral_radius = toteig #average
standard_deviation = sqrt(toteigsq - toteig*toteig) #standard deviation


################
#Writing results
################

results_bis = "spectral_radius.txt"
R = mat.or.vec(1,4)
R = rbind(intw, real_spectral_radius, average_spectral_radius, standard_deviation)
write.matrix(t(R), file=results_bis, sep=" ")

