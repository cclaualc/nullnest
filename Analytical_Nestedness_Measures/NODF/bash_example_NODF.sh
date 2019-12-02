
#We start compiling the program using the compiler gfortran
gfortran -Og -o "NODF_analytic.exe" NODF_analytic.f90

#We run the program for the example files indexed i=0
./NODF_analytic.exe 0 
