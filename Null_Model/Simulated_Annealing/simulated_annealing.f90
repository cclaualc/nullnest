
      module global
            integer :: N_R, N_C !dimension of the bipartite matrix
            real*8, allocatable :: nr_d(:), nc_d(:) !vectors containing the degree sequences of rows and columns
      end module global


!*********************************************************************************************
!*******************************PROGRAM: SIMULATED ANNEALING**********************************

      !This program constructs a null model for bipartite networks, based on a maximum-entropy ensemble that preserves degree sequences on average.
      !It finds the Lagrange multipliers that maximize the likelihood of finding the real degree sequences, then writes the probability matrix of interactions.
      !The optimization method is a simulated annealing algorithm.

      !Code by Clàudia Payrató Borràs [1], partly based on the simulated annealing code by Goffe et al. [2]
      !For any questions on the code, you may contact me by mail: claudiapb13@gmail.com

      !When using this program please cite: 
      ![1] "Breaking the Spell of Nestedness: "The Entropic Origin of Nestedness in Mutualistic Systems", C. Payrató-Borràs et al. 
      !Phys. Rev. X 9 (2019)
      ![2] "Global Optimization of Statistical Functions with Simulated Annealing," Goffe et al. Journal of Econometrics, vol. 60 (1994)
      ![3] "Generation of Gaussian distributed random numbers by using a numerical inversion method, R. Toral, A. Chakrabarti, Computer
      !Physics Communications 74, 327 (1993)

!*********************************************************************************************

      PROGRAM SA
      use global
      implicit none
      character (len=25) arg(2), rstr, generaldoc, matrixdoc, mranddoc, lagrangedoc, maxdoc
      integer :: I, J, K, inew, irep, idoc, iloop
      integer :: N, N_T, N_S, ndown, nup, Nredim
      integer, allocatable :: M(:,:), nacp(:), point(:), cont(:)
      real*8, allocatable :: x(:), M_test(:,:), c(:), xopt(:), xp(:), fstar(:)
      real*8, allocatable :: nr_d_test(:), nc_d_test(:), lb(:), ub(:), vm(:)
      real*8 f, fp, fopt, T, p, ran, rbuff, ratio, eps, rt, deg, delta_ran, xtest, prod
      external FCN
      real*8 dran_u
      integer i_dran
      logical quit, loop
      quit = .false.
      irep = 5 !Number of convergent iterations requested


!*******************************************************
!******DEFINITION AND INITIALIZATION of VARIABLES******
!*******************************************************

!****External arguments (to be entered by terminal)
      !Network index 'idoc', defines the real network upon which the null model is constructed
      call getarg(1, arg(1))
      read(arg(1),*) idoc
      !Index of the iteration index 'iloop', to be changed if one aims to perform independent runs of the simulated annealing
      call getarg(2, arg(2))
      read (arg(2),*) iloop

!****Opening files 
      write (rstr, '(I3.1)') idoc
      generaldoc='general'//trim(adjustl(rstr))//'.txt' !Input file containing dimension
      matrixdoc='matrix'//trim(adjustl(rstr))//'ord.txt' !Input file containing ordered real matrix
      mranddoc='matrix'//trim(adjustl(rstr))//'rand.txt' !Output file containing probability matrix
      maxdoc='max_loglikelihood'//trim(adjustl(rstr))//'.txt' !Output file containing final value of the optimized log-likelihood
      lagrangedoc='lagrange_multipliers'//trim(adjustl(rstr))//'.txt'!Output file containing the final values of the Lagrange Multipliers
      open (10,file=generaldoc,status='old')
      open (11, file=matrixdoc, status='old')

!*****Reading files and allocating variables
      read(10,*) N_R, N_C !Number of rows, number of columns
      N = N_R+N_C !Number of unknown variables and thus equations that we need to solve
      allocate (x(N), M(N_R, N_C), M_test(N_R, N_C), nr_d(N_R), nc_d(N_C), nr_d_test(N_R), nc_d_test(N_C))
      allocate (lb(N), ub(N), xopt(N), xp(N), fstar(irep), point(N))
      
!*****Definition and initialization of the variables
      rt = 0.85 !Reduction factor of the temperature
      f = -1.d100 !f is the trial value of the function to be optimized
      fopt = -1.d100 !fopt is the proposed global maximum
      T = 1000.d0 !Initial temperature 
      x(1:N) = 1.d0 !Initial values of the vector containing the Lagrange multipliers (x)
      N_T = 200 !Number of updates per fixed temperature
      N_S = 10 !Number of updates per fixed exploring distance 'vm'
      lb(1:N) = 1.d-20 !lower bound for the Lagrange multipliers
      ub(1:N) = 1.d20 !upper bound for the Lagrange multipliers
      eps = 10.d-6 !Convergence tolerance
      fstar(1:irep) = 1.d20 !Vector to be used to check the convergence criteria
      deg = 0 !Degeneracy of degrees


!****Reading the real matrix M and calculating the degree sequences 
      !Empirical degree sequences of rows in vector nr_d
      do I=1,N_R
            nr_d(I) = 0.d0
            !Reading M
            read(11,*) (M(I,J), J=1,N_C)
            do J=1,N_C
                  nr_d(I) = M(I,J) + nr_d(I)
            enddo
      enddo
      !Empirical degree sequences of columns in vector nc_d
      do J=1,N_C
            nc_d(J) = 0.d0
            do I=1,N_R,1
                  nc_d(J) = M(I,J) + nc_d(J)
            enddo
      enddo

      close(10)
      close(11)


!****Redimensionalization of the system to get ride of equivalent equations 
      !This reduces the dimensionality of the search space and enhances the optimization.
      !*REMINDER*: you need to input a SORTED matrix by decreasing degree.
      K = 1
      point(1) = 1
      I = 1
      do while (I .lt. N_R)
            J = I + 1
            loop = .true.
            do while ((loop .eqv. .true.).and.(J .le. N_R))
                  if (nr_d(I) .eq. nr_d(J)) then
                        point(J) = K
                        deg = deg + 1
                        J = J + 1
                        loop = .true.
                  else
                        K = K + 1
                        point(J) = K
                        loop = .false.
                  endif
            enddo
            I = J
      enddo
      I = 1
      K = K + 1
      point(N_R+1) = K
      do while (I .lt. N_C) 
            J = I+1
            loop = .true.
            do while ((loop .eqv. .true.) .and. (J .le. N_C))
                  if (nc_d(I) .eq. nc_d(J)) then
                        point(J + N_R) = K
                        deg = deg + 1
                        J = J + 1
                        loop = .true.
                  else
                        K = K + 1
                        point(J + N_R) = K
                        loop = .false.
                  endif
            enddo
            I = J
      enddo
      !Nredim = number of nodes - number of repeated degrees
      Nredim = N - deg
      !We allocate the necessary vectors
      allocate(cont(Nredim), nacp(Nredim), vm(Nredim))
      vm(1:Nredim) = 1.d0
      N_S = N_S*Nredim

!*****Initialize the random number generator from external random number generator dranxor
      call dran_ini(12345+iloop)


!**************************************************
!******OPTIMIZATION using SIMULATED ANNEALING******
!*************************************************

!     Here we find the Lagrange Multipliers that optimize the log-likelihood of finding degree sequences in the ensemble
!     The simulated annealing algorithm perfoms a pseudo-random search along the space of parameters
!     It accepts points following the Metropolis criteria. This allows for downhill moves and thus to scape from local maxima
!     The present implementation is an adaptation to our particular problem of the code provided by Goffe, Ferrier and Rogers [2],
!     which on turn is based on the simulated annealing algorithm proposed by Corana et al. (see Documentation for more details)

!*****First loop over T values
      !In this loop we progressively reduce the value of the temperature, T.
      !This procedure progressively reduces the probability of accepting downhill moves. 
      !The total number of trial points explored at each temperature is: N_S*N_T (where N_S is proportional to Nredim).
      !The loop finishes when the algorithm has converged to a global maximum by fulfilling the convergence criteria
      do while (quit .eqv. .false.)

!*****Second loop with varying 'vm'. In this loop we adjust the value of 'vm' (exploring distance) in order to ensure an exhaustive exploration of the search space.
      do I=1, N_T, 1

            nacp(1:Nredim) = 0 !number of accepted points
            cont(1:Nredim) = 0 !total number of proposed points

            !****Third loop with fixed 'vm'
            do J=1, N_S,1

                  !Pick randomly one of the elements of xp, which we will change and save the index as 'inew'
                  xp(1:N) = x(1:N)
                  inew = i_dran(Nredim)

                  !Produce a random distance 'delta_ran'
                  delta_ran = 2.d0*dran_u() - 1.d0
                  cont(inew) = cont(inew) + 1

                  !Change all equivalent nodes to 'inew' (having the same degree, on the same guild)
                  do K=1, N,1
                         if (point(K) .eq. inew) then
                              xp(K) = x(K) + vm(inew)*delta_ran
                              xtest = xp(K)
                         endif
                  enddo

                  !Check if xtest is inside the bounds, otherwise find another one which falls inside.
                  do while ((xtest < lb(inew)) .or. (xtest > ub(inew))) 
                        delta_ran = 2.d0*dran_u() - 1.d0
                        do K=1, N,1
                              if (point(K) .eq. inew) then
                                    xp(K) = x(K) + vm(inew)*delta_ran
                                    xtest = xp(K)
                              endif
                        enddo
                  enddo

                  !Evaluation of the function to optimize at the trial point xp (here, the log-likelihoog). The subroutine returns fp.
                  call FCN(N,xp,fp)

                  !****Acceptance/rejection of the trial point
                  !If the function 'fp' is greater than the previous value 'f', accept xp
                  if (fp .gt. f) then
                        x(1:N) = xp(1:N)
                        f = fp
                        nacp(inew) = nacp(inew) + 1 !new acceptance
                        !If 'fp' is larger than the highest maximum recorded up to now 'fopt', set 'fp' as the new global maximum
                        if (fp .gt. fopt) then
                              fopt = fp
                              xopt(1:N) = x(1:N) !We save the coordinates of fopt
                        endif
                  !If the function value 'fp' is lower than 'f', accept/reject following the Metropolis criteria
                  else
                        p = (fp - f)/T 
                        ran = dlog(dran_u())
                        if (ran .le. p) then 
                              x(1:N) = xp(1:N)
                              f = fp
                              nacp(inew) = nacp(inew) + 1 !new acceptance
                        endif
                  endif

            enddo !End of the loop with fixed vm

            !Change of vm to accept aproximately half of the trials. This warrants that the search space is properly explored.
            do K = 1, Nredim, 1
                  if (cont(K) > 0) then
                        ratio = dble(nacp(K))/dble(cont(K))
                        if (ratio .gt. 0.6d0) then
                              vm(K) = vm(K)*(1.d0 + (ratio - 0.6d0)/0.8d0)
                        elseif (ratio .lt. 0.4d0) then
                              vm(K) = vm(K)/(1.d0 + (0.4d0 - ratio)/0.8d0)
                        endif
                        if (vm(K) .gt. (ub(K)-lb(K))) then
                              vm(K) = (ub(K) - lb(K))/10.d0
                        elseif (vm(K) .lt. 1.d-4) then
                              vm(K) = 1.d0
                        endif
                  endif
            enddo


      enddo !End of the loop with fixed temperature
  
      !Check if the convergence criteria are fulfilled: 
      !At least the last 'irep' runs must have converged to a function value differing less than 'eps'
      fstar(1) = f
      if (abs(fopt - fstar(1)) <= eps) quit = .true.
      do K = 2, irep,1
            if (abs(f - fstar(K)) > eps) quit = .false.
      enddo
      !Reduce the temperature and set the starting point x equal to the previous optimum xopt, so the next iteration will start in the most promising region
      T = rt*T
      f = fopt
      x(1:N) = xopt(1:N)

      !Update the vector 'fstar' containing the resulting fopt of the last iterations
      do K=irep, 2, -1
            fstar(K) = fstar(K-1) 
      enddo

      !*If you'd like to see the evolution of the search by printing the current temperature and global maximum on the terminal, un-comment the following command*
      !write(*,*) "T=",T,"Global maximum found up to now=",fopt

      enddo !End of the general do-while loop


!**************************************************
!******TEST & WRITING OF FINAL RESULTS******
!*************************************************

!****Construct the matrix containing the probability of interaction among species, 'M_test'
      do I=1, N_R, 1
            do J=1, N_C, 1
                  prod = xopt(I)*xopt(J+N_R)
                  M_test(I,J) = prod/(1.d0+prod)
            enddo
      enddo

!*****Test the results by checking that the degree sequences are kept equal to the real ones
      do I=1,N_R
            nr_d_test(I) = 0.d0
            do J=1,N_C
                  rbuff = nr_d_test(I)
                  nr_d_test(I) = M_test(I,J) + rbuff
            enddo
      enddo
      do J=1,N_C
            nc_d_test(J) = 0.d0
            do I=1,N_R
                  nc_d_test(J) = M_test(I,J) + nc_d_test(J)
            enddo
      enddo
      !*To perform a comparison among the real degree sequences and the average ones over the ensemble, uncomment the following lines*
      !write(*,*) "Differences between real degrees and average degrees, for rows="
      !do I=1, N_R, 1
      !      write(*,*) abs(nr_d(I) - nr_d_test(I))
      !enddo
      !write(*,*) "Differences between real degrees and average degrees, for columns="
      !do I=1, N_C, 1
      !      write(*,*) abs(nc_d(I) - nc_d_test(I))
      !enddo


!*****Writing the results on external files

      !Global maximum for the log-likelihood
      open(18, file=maxdoc, status='unknown') !Comment this line to avoid overwriting previous results
      !*In order to save the present result without over-writing the result of previous runs, comment the previous line and uncomment the next line*
      !open(18, file=maxdoc, status='old',position='append') !Attention: the file we're opening needs to exist already
      write(18,*) "Run index= ",iloop, "         Global maximum of log-likelihood=", fopt
      close(18)

      !Probability matrix of interactions
      open(14, file=mranddoc, status='unknown') !Comment this line to avoid overwriting previous results
      !*In order to save the present result without over-writing the result of previous runs, comment the previous line and uncomment the next line*
      !open(14, file=mranddoc, status='old',position='append') !Attention: the file we're opening needs to exist already
      do I=1, N_R
            write(14,*) (M_test(I,J), J=1,N_C) 
      enddo
      close(14)

      !Optimal Lagrange multipliers
      open(15, file=lagrangedoc, status='unknown') !Comment this line to avoid overwriting previous results
      !*In order to save the present result without over-writing the result of previous runs, comment the previous line and uncomment the next line*
      !open(15, file=lagrangedoc, status='old',position='append') !Attention: the file we're opening needs to exist already
      write(15,*) "Lagrange multipliers for rows="
      do I=1, N_R, 1
            write(15,*) xopt
      enddo
      write(15,*) "Lagrange multipliers for columns="
      do I=N_R+1, N_R+N_C, 1
            write(15,*) xopt
      enddo
      close(15)

      end program SA
!****************************************************
!****************************************************



!**************RANDOM NUMBERS GENERATOR*****
!******External random number generator 'dranxor.f90', by Raul Toral (see reference [3] above)************
!******Please when using this code ackownledge the authors by citing [3]
      INCLUDE "dranxor.f90"
!******



!****************************************************
!*****************SUBROUTINE FCN*********************
!****************************************************

!*****This subroutine calculates the loglikelihood that will be maximised
!*****As inputs, it takes the number of variables, N, and the trial vector of Lagrange multipliers (x) 
!*****As output, it provides the value of the log-likelihood, f
!*****Uses as well the global variables N_R, N_C, nr_d and nc_d

      subroutine FCN(N,x,f)
            use global
            implicit none
            integer :: N, K, J
            real*8 :: x(N)
            real*8 :: f, L, lbuff
            L = 0.d0
            !Sum over rows' contributions to the log-likelihood
            do K=1, N_R, 1
                  lbuff = L
                  L = lbuff + nr_d(K)*dlog(x(K))
                  do J=N_R+1,N, 1
                        lbuff = L
                        L = L - dlog(1.d0+x(K)*x(J))
                  enddo
            enddo
            !Sum over columns' contributions to the log-likelihoods
            do K=N_R+1, N, 1
                  L = L + nc_d(K-N_R)*dlog(x(K))
            enddo
            f = L
            return
      end subroutine fcn
!****************************************************

