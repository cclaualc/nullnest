library(MASS)
library('bipartite')


##########################PROGRAM Temperature########################
#This program calculates a standarized version of the nestedness metrics Temperature using the R package 'bipartite'
#It produces the following information:
      #i) The temperature of a real network provided by the user
      #ii) The average temperature over a sampling of a null ensemble
      #iii) The standard deviation of the distribution of the temperature metrics over the ensemble
##The program takes as input the ordered real bi-adjacency matrix as well as the probability matrix over the null ensemble
#Code made by Clàudia Payrató Borràs. For any questions on the code, you may contact me by mail: claudiapb13@gmail.com
#When using this program, please cite: 
#[1] "Lost in Nestedness? Assessing the performance of different metrics" Payrató-Borràs et al. (2019)
#[2] "Introducing the bipartite package: analysing ecological networks." Dormann et al.Interaction 1.0.2413793 (2008).


##Input arguments
args <- commandArgs(TRUE)
#Index 'intw', which identifies the network to study
intw <- as.integer(args[1])
#Size of the sampling: number of null networks to generate. We recommend a minimum size of 10³ networks
N_sampl <- as.integer(args[2])


########################################################
########Calculation of the REAL Temperature##########
########################################################

#Opening the file containing the probability matrix
file_real = paste("matrix",intw,"ord.txt",sep="")

#Reading the real bi-adjacency matrix m
m = read.table(file_real)
m = as.matrix(m)

#Calculation of the real temperature
tm <- nested(m, method = "binmatnest2", rescale=FALSE, normalised=FALSE)

#Standarization of the temperature, so that larger temperature means larger nestedness
temp = 100 - tm
real_temperature = temp



########################################################
#####SAMPLING of the Temperature over the ensemble######
########################################################

#Opening the file containing the probability matrix
file = paste("matrix",intw,"rand.txt", sep="")

#Reading the probability matrix of interactions, p
p = read.table(file)

#Initializing vectors and parameters
set.seed(12345) #initializing random number generator
tottemp = 0.0 #sum of the temperature of each null network
tottempsq = 0.0 #sum of the squared temperature of each null network

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

#Calculation of the standarized temperature
tm <- nested(b, method = "binmatnest2", rescale=FALSE, normalised=FALSE)
temp = 100 - tm #We invert it, so that larger temperature means more nested

#Calculation of the average and squared average of the temperature
tottemp = tottemp + temp/N_sampl
tottempsq = tottempsq + temp*temp/N_sampl
iter = iter + 1

} #end of the sampling iteration
 
#Final calculations
average_temperature = tottemp #average
standard_deviation = sqrt(tottempsq - tottemp*tottemp) #standard deviation


################
#Writing results
################

results_bis = "temperature.txt"
R = mat.or.vec(1,4)
R = rbind(intw,real_temperature,average_temperature,standard_deviation)
write.matrix(t(R), file=results_bis, sep=" ")

