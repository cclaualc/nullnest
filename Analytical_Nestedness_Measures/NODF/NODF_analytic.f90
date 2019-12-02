      module global
            integer :: N_R, N_C !Global parameters: number of rows, number of columns
      end module global



!*********************************************************************************************
!***********************************PROGRAM: NODF_analytic********************************

!     This program calculates, using the metrics NODF:
              !i) the real value of nestedness of an empirical network
              !ii) the first two moments of the metric, in a maximum entropy ensemble, by using the analytical expressions derived in [1]
!     It takes as input the ordered real bi-adjacency matrix as well as the probability matrix over the null ensemble
!     It returns as output the real measure of NODF, the average expected value in the ensemble and the standard deviation in the ensemble

!     Code by Clàudia Payrató Borràs. For any questions on the code, you may contact me by mail: claudiapb13@gmail.com
!     When using this program, please cite: 
!     [1] "Breaking the Spell of Nestedness: The Entropic Origin of Nestedness in Mutualistic Systems", C. Payrató-Borràs et al. 
!     Phys. Rev. X 9 (2019)

!*********************************************************************************************


      program NODF
      use global
      implicit none
      character(len=25) arg(1), generaldoc, matrixdoc, mranddoc, resultsdoc, rstr
      integer :: N, I, J, K, Nprow, Npcol, idoc, iter
      integer, allocatable :: M(:,:)
      real*8, allocatable :: Mrand(:,:), Npaired_col(:), Npaired_row(:)
      real*8, allocatable :: nr_d(:), nc_d(:)
      real*8, allocatable :: nr_drand(:), nc_drand(:)
      real*8 :: nest_real, suma, nest_rand, sigma, sigma_ij, derivative, z, rel_d


!*******************************************************
!******DEFINITION AND INITIALIZATION of VARIABLES******
!*******************************************************

!****External arguments (to be entered by terminal)
      !Network index 'idoc', defines the network for which we calculate the NODF measures
      call getarg(1, arg(1))
      read(arg(1),*) idoc

!****Opening files 
      write (rstr, '(I3.1)') idoc
      generaldoc='general'//trim(adjustl(rstr))//'.txt'
      matrixdoc='matrix'//trim(adjustl(rstr))//'ord.txt'
      mranddoc='matrix'//trim(adjustl(rstr))//'rand.txt'
      resultsdoc='analytic_NODF.txt'
      open (10,file=generaldoc,status='old')
      open (11, file=matrixdoc, status='old')
      open(12, file=mranddoc, status='old')

!*****Reading files and allocating variables
      read(10,*) N_R, N_C !Number of rows, number of columns
      N = N_R+N_C !Total number of scpecies
      Nprow = (N_R*(N_R-1))/2 !Number of pairs of rows
      Npcol = (N_C*(N_C-1))/2 !Number of pairs of columns
      allocate (M(N_R, N_C), Mrand(N_R, N_C), nr_d(N_R), nc_d(N_C))
      allocate (nr_drand(N_R), nc_drand(N_C),Npaired_row(Nprow), Npaired_col(Npcol))

!*****Calculating the degree sequences and reading the matrices
      do I=1,N_R
            nr_d(I) = 0.d0
            nr_drand(I) = 0.d0
            read(11,*) (M(I,J), J=1,N_C) !Reading the real bi-adjacency matrix
            read(12,*) (Mrand(I,J), J=1,N_C) !Reading the probability bi-adjacency matrix
            !Degree sequences of rows: nr_d (real matrix) and nr_drand (average in the ensemble)
            do J=1,N_C
                  nr_d(I) = dble(M(I,J)) + nr_d(I)
                  nr_drand(I) = Mrand(I,J) + nr_drand(I)
            enddo
      enddo

      do J=1,N_C
            nc_d(J) = 0.d0
            nc_drand(J) = 0.d0
            !Degree sequences of columns: nc_d (real matrix) and nc_drand (average in the ensemble)
            do I=1,N_R
                  nc_d(J) = dble(M(I,J)) + nc_d(J)
                  nc_drand(J) = Mrand(I,J) + nc_drand(J)
            enddo
      enddo


!**************************************************
!*****REAL NETWORK NESTEDNESS CALCULATION**********
!**************************************************

      !Call of the subroutine that calculates the unnormalised nestedness among columns
      call nested_col(dble(M), nc_d, Npcol, Npaired_col)

      !Call of the subroutine that calculates the unnormalised nestedness amongg rows
      call nested_row(dble(M), nr_d, Nprow, Npaired_row)

      !Calculation of the real global nestedness
      nest_real = 0.d0
      do I=1, Npcol, 1
            suma = nest_real
            nest_real = suma + Npaired_col(I)/dble(Npcol+Nprow)
      enddo
      do I=1, Nprow, 1
            suma = nest_real
            nest_real = suma + Npaired_row(I)/dble(Npcol+Nprow)
      enddo


!**************************************************
!*****RANDOMISED NETWORK NESTEDNESS CALCULATION****
!**************************************************

     !Call of the subroutine that calculates the unnormalised nestedness among columns
      call nested_col(Mrand, nc_drand, Npcol, Npaired_col)

     !Call of the subroutine that calculates the unnormalised nestedness among rows
      call nested_row(Mrand, nr_drand, Nprow, Npaired_row)

     !Calculation of the average global nestedness
      nest_rand = 0.d0
      do I=1, Npcol, 1
            suma = nest_rand
            nest_rand = suma + Npaired_col(I)/dble(Npcol+Nprow)
      enddo
      do I=1, Nprow, 1
            suma = nest_rand
            nest_rand = suma + Npaired_row(I)/dble(Npcol+Nprow)
      enddo

     !Calculation of the standard deviation 
      sigma = 0
      do I=1, N_R
            do J=1, N_C
                  sigma_ij = Mrand(I,J)*(1.d0-Mrand(I,J))
                  call derivative_sigma(Mrand, I, J, nr_d, nc_d, Nprow, Npcol, derivative)
                  sigma = sigma + derivative*derivative*sigma_ij
            enddo
      enddo
      sigma = dsqrt(sigma)


!**************************************************
!******WRITING FINAL RESULTS******
!*************************************************

      open(13, file=resultsdoc, status='unknown')
      write(13,*) "Index of the network=",idoc,"Real NODF=",nest_real,"Average NODF=",nest_rand,"Standard deviation=",sigma

      close(10)
      close(11)
      close(12)

      end program NODF  
!************************************************************************



!****************************************************
!**************SUBROUTINE nested_col****************
!****************************************************
!****This subroutine calculates the unnormalised contribution to NODF of all pairs of columns**** 

      subroutine nested_col(M,MT,Npcol,Npaired_col)
      use global
      implicit none
      integer :: PO, I, J, K, L, Npcol
      real*8 :: Npaired_col(Npcol), M(N_R,N_C), MT(N_C), buff, suma 
      K=1

      do I=1, N_C-1, 1
            do J=I+1,N_C,1
                  !***Decreasing fill condition (DF)
                  if (MT(I).gt. MT(J)) then
                        !DF=100, we take into account the contribution
                        suma = 0.d0
                        do L=1, N_R, 1
                              buff = suma
                              suma = buff + dble(M(L,I))*dble(M(L,J))
                        enddo
                        Npaired_col(K) = 100.d0*dble(suma)/dble(MT(J))
                  elseif (MT(I) .le. MT(J)) then
                        !DF=0, we neglect the contributions
                        Npaired_col(K)=0
                  endif
                  K=K+1
            enddo 
      enddo
      end


!****************************************************
!**************SUBROUTINE nested_row**************
!****************************************************
!****This subroutine calculates the unnormalised contribution to NODF of all pairs of rows**** 

      subroutine nested_row(M,MT,Nprow,Npaired_row)
      use global
      implicit none
      integer :: PO, I, J, K, L, Nprow
      real*8 :: Npaired_row(Nprow), M(N_R,N_C), MT(N_R), buff, suma 
      K=1
      do I=1, N_R-1, 1
            do J=I+1,N_R,1
                  !***Decreasing fill condition (DF)
                  if (MT(I).gt. MT(J)) then
                        !DF=100
                        suma = 0
                        do L=1, N_C, 1
                              buff = suma
                              suma = buff + dble(M(I,L))*dble(M(J,L))
                        enddo
                        Npaired_row(K) = 100.d0*dble(suma)/dble(MT(J))
                  elseif (MT(I) .le. MT(J)) then
                        !DF=0
                        Npaired_row(K)=0
                  endif
                  K=K+1
            enddo 
      enddo
      end


!****************************************************
!**************SUBROUTINE derivative_sigma***********
!****************************************************
!****This subroutine calculates the derivative of NODF with respect to a matrix element m_{ab}, using the expression obtained in [1]**** 

      subroutine derivative_sigma(M,a,b, MT_row, MT_col, Nprow, Npcol, derivative)
      use global
      implicit none
      integer:: a, b, I, J, K, L, Nprow, Npcol
      real*8 :: M(N_R, N_C), MT_row(N_R), MT_col(N_C)
      real*8 :: derivative
      derivative = 0.d0

     !***Derivative of the ROW term
      if (a .lt. N_R) then
            do J = a+1,N_R,1
                  if (MT_row(a) .gt. MT_row(J)) then !DF = 100
                        derivative = derivative + 100.d0*M(J,b)/MT_row(J)
                  endif
            enddo
      endif
      do I = 1, a-1, 1
            if (MT_row(I) .gt. MT_row(a)) then !DF=100
                  derivative = derivative + 100.d0*M(I,b)/MT_row(a)
            endif 
      enddo
      do I=1, a-1,1
            do K =1, N_C, 1
                  if (MT_row(I) .gt. MT_row(a)) then !DF=100
                        derivative = derivative - 100.d0*M(I,K)*M(a,K)/(MT_row(a)*MT_row(a))
                  endif
            enddo
      enddo

     !***Derivative of the COLUMN term
      if (b .lt. N_C) then 
            do L = b+1,N_C,1
                  if (MT_col(b) .gt. MT_col(L)) then !DF=100
                        derivative = derivative +100.d0*M(a,L)/MT_col(L)
                  endif
            enddo
      endif
      do K=1, b-1,1
            if (MT_col(K) .gt. MT_col(b)) then !DF=100
                  derivative = derivative + 100.d0*M(a,K)/MT_col(b)
            endif
      enddo
      do K=1,b-1,1
            do I=1, N_R,1
                  if (MT_col(K) .gt. MT_col(b)) then !DF=100
                       derivative = derivative - 100.d0*M(I,K)*M(I,b)/(MT_col(b)*MT_col(b))
                  endif
            enddo
      enddo

      !Normalization of the result
      derivative = derivative/dble(Nprow+Npcol)
      return

      end
