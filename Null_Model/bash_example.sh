
#We start compiling the program using the compiler gfortran
gfortran -Og -o "simulated_annealing.exe" simulated_annealing.f90

#We run the program for the example files indexed i=0 and a first run indexed iloop=1
./simulated_annealing.exe 0 1
