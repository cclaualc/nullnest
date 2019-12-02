
#We start compiling the program using the compiler gfortran
gfortran -Og -o "NIR.exe" NIR.f90

#We run the program for the example files indexed i=0 and generate a sample size of 10² null networks
./NIR.exe 0 100
