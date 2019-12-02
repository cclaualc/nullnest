library('bipartite')
library(MASS)

##########################PROGRAM Discrepancy########################
#This program calculates the nestedness metrics Discrepancy using the R package 'bipartite'
#It produces the following information:
      #i) The discrepancy of a real network provided by the user
      #ii) The average discrepancy over a sampling of a null ensemble
      #iii) The standard deviation of the distribution of the discrepancy over the null ensemble
##The program takes as input the ordered real bi-adjacency matrix as well as the probability matrix over the null ensemble
#Code made by Clàudia Payrató Borràs. For any questions on the code, you may contact me by mail: claudiapb13@gmail.com
#When using this program, please cite: 
#[1] "Lost in Nestedness? Assessing the performance of different metrics" Payrató-Borràs et al., pre-rpint (2019).
#[2] "Introducing the bipartite package: analysing ecological networks." Dormann et al., Interaction 1.0.2413793 (2008).


##Input arguments by terminal
args <- commandArgs(TRUE)
#Index 'intw', which identifies the network to study
intw <- as.integer(args[1])
#Size of the sampling 'N_sampl': number of null networks to generate. We recommend a minimum size of 10³ networks
N_sampl <- as.integer(args[2])


########################################################
########Calculation of the REAL discrepancy#############
########################################################

#Opening the file containing the probability matrix
file_real = paste("matrix",intw,"ord.txt",sep="") 

#Reading the real bi-adjacency matrix m
m = read.table(file_real)
m = as.matrix(m)

#Calculation of the discrepancy of m
a <- nested(m, method = "discrepancy", rescale=FALSE, normalised=FALSE)

#Calculation of the fill (total number of links)
conn = networklevel(m, index="connectance", level="both", weighted=FALSE)
fill = conn*ncol(m)*nrow(m)

#Standarization of the discrepancy by dividing by the fill, inverting it and rescaling it between 0 and 100
real_discrepancy = 100 - 100*a/fill 


###############################################################
########SAMPLING of the Discrepancy over the ensemble##########
###############################################################

#Opening the file containing the probability matrix
file_rand = paste("matrix",intw,"rand.txt",sep="") 

#Reading the probability matrix of interactions, p
p = read.table(file_rand)

#Initializing vectors and parameters
set.seed(123456) #initializing random number generator
totdis = 0.0 #sum of the discrepancy of each null network
totdissq = 0.0 #sum of the squared discrepancy of each null network

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

#Calculation of the non-standarized discrepancy
a <- nested(b, method = "discrepancy", rescale=FALSE, normalised=FALSE)

#Calculation of the fill (total number of links)
conn = networklevel(b, index="connectance", level="both", weighted=FALSE)
fill = conn*ncol(b)*nrow(b)

#Calculation of the standarized version of the discrepancy (normalized by the density of links and inverted)
dis = 100 - 100*a/fill #We invert it, so that larger discrepancy means more nested

#Calculation of the average and the squared average of the discrepancy
totdis = totdis + abs(dis/N_sampl)
totdissq = totdissq + abs(dis*dis/N_sampl)
iter = iter + 1

} #end of the sampling iteration

#Final calculations
average_discrepancy = totdis
standard_deviation = sqrt(totdissq - totdis*totdis)


################
#Writing results
################
results_bis = "discrepancy.txt"
R = mat.or.vec(1,4)
R = rbind(intw,real_discrepancy,average_discrepancy,standard_deviation)
write.matrix(t(R), file=results_bis, sep=" ")

