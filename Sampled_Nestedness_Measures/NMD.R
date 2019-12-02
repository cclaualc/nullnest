library('bipartite')
library(MASS)


##########################PROGRAM MANHATTAN DISTANCE########################
#This program calculates a standarized version of the nestedness index based on Manhattan distance (NMD from now on), using the R package 'bipartite'
#It produces the following information:
      #i) The NMD of a real network provided by the user
      #ii) The average NMD over a sampling of a null ensemble
      #iii) The standard deviation of the distribution of the NMD metrics over the ensemble
##The program takes as input the ordered real bi-adjacency matrix as well as the probability matrix over the null ensemble
#Code made by Clàudia Payrató Borràs. For any questions on the code, you may contact me by mail: claudiapb13@gmail.com
#When using this program, please cite: 
#[1] "Breaking the Spell of Nestedness: "The Entropic Origin of Nestedness in Mutualistic Systems", C. Payrató-Borràs et al. Phys. Rev. X 9 (2019)
#[2] "Introducing the bipartite package: analysing ecological networks." Dormann et al., Interaction 1.0.2413793 (2008).


##Input arguments by terminal
args <- commandArgs(TRUE)
#Index 'intw', which identifies the network to study
intw <- as.integer(args[1])
#Size of the sampling 'N_sampl': number of null networks to generate. We recommend a minimum size of 10³ networks
N_sampl <- as.integer(args[2])
#Number of null configuration 'N_reps' produced in the calculation of NMD, randomized keeping constant the size and fill 
N_reps <- as.integer(args[3])


########################################################
########Calculation of the REAL NMD#############
########################################################

#Opening the file containing the probability matrix
file_real = paste("matrix",intw,"ord.txt",sep="") 

#Reading the real bi-adjacency matrix 'm'
m = read.table(file_real)
m = as.matrix(m)

#Calculation of the NMD of 'm'
nmd <- nestedness.corso(m, weighted=FALSE, reps=N_reps) #The parameter 'reps' represents the number of random matrices drawn to provide the significant measure. 

#Standarization of the NMD by dividing by the fill, inverting it and rescaling it between 0 and 100
real_nmd = 100.0*(1.0 - nmd)


###############################################################
########SAMPLING of the NMD over the ensemble##########
###############################################################

#Opening the file containing the probability matrix
file_rand = paste("matrix",intw,"rand.txt",sep="") 

#Reading the probability matrix of interactions, p
p = read.table(file_rand)

#Initializing vectors and parameters
set.seed(123456) #initializing random number generator
totnmd = 0.0 #sum of the NMD of each null network
totnmdsq = 0.0 #sum of the squared NMD of each null network

#Start of the iteration loop to produce the sample
iter = 0
while (iter < N_sampl){

#Generation of a new null network b, using the matrix probability p
b = runif(nrow(p)*ncol(p))
b = matrix(b,nrow(p),ncol(p))
aux = b>p
b[aux]=0
b[!aux] = 1
b = as.matrix(b)

#Calculation of the non-standarized NMD
a <- nestedness.corso(b, weighted=FALSE, reps=N_reps) #The parameter 'reps' represents the number of random matrices drawn to provide the significant measure.

#Standarization of the NMD
nmd = 100.0*(1.0 - a) #We invert it, so that larger NMD means more nested

#Calculation of the average and the squared average of the NMD
totnmd = totnmd + abs(nmd/N_sampl)
totnmdsq = totnmdsq + abs(nmd*nmd/N_sampl)
iter = iter + 1
} #end of the sampling iteration

#Final calculations
average_nmd = totnmd #average
standard_deviation = sqrt(totnmdsq - totnmd*totnmd) #standard deviation


################
#Writing results
################

results_bis = "NMD.txt"
R = mat.or.vec(1,4)
R = rbind(intw,real_nmd,average_nmd,standard_deviation)
write.matrix(t(R), file=results_bis, sep=" ")

