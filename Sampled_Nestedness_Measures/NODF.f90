      module global
            integer :: N_R, N_C !Global variables: number of rows, number of columns
      end module global


!*********************************************************************************************
!***********************************PROGRAM: NODF********************************

!     This program calculates, for the nestedness metrics NODF:
              !i) the real value of nestedness of an empirical network
              !ii) the first two moments of the metric, in a maximum entropy ensemble, computed by numerically sampling the ensemble
!     It takes as input the ordered real bi-adjacency matrix as well as the probability matrix over the null ensemble
!     It returns as output the real measure of NODF, the average expected value in the ensemble and the standard deviation in the ensemble

!     Code by Clàudia Payrató Borràs. For any questions on the code, you may contact me by mail: claudiapb13@gmail.com
!     When using this program, please cite: 
!     [1] "Breaking the Spell of Nestedness: The Entropic Origin of Nestedness in Mutualistic Systems", C. Payrató-Borràs et al. 
!     Phys. Rev. X 9 (2019)
!     [2] "Generation of Gaussian distributed random numbers by using a numerical inversion method, R. Toral, A. Chakrabarti, Computer
      !Physics Communications 74, 327 (1993)

!*********************************************************************************************


      program NODF
      use global
      implicit none
      character(len=25) arg(2), generaldoc, matrixdoc, mranddoc, resultsdoc, rstr
      integer :: N, I, J, K, Nprow, Npcol, idoc, iter, Nsampl, izeros_rows, izeros_col, N_Cclean, N_Rclean
      integer, allocatable :: M(:,:), Mrand(:,:), Mrand0(:,:), Mclean(:,:)
      real*8, allocatable :: Prand(:,:), Npaired_col(:), Npaired_row(:), nodf_vec(:,:)
      integer, allocatable :: nr_d(:), nc_d(:), nr_drand(:), nc_drand(:), nr_drandclean(:), nc_drandclean(:)
      real*8 :: rbuff, rbuffrand, nest_real, suma, nest_rand, sigma, sigma_ij, derivative, z, rel_d
      real*8 :: dran_u, u, suma2


!*******************************************************
!******DEFINITION AND INITIALIZATION of VARIABLES******
!*******************************************************

!****External arguments (to be entered by terminal)
      !Network index 'idoc', defines the real network upon which the null model is constructed
      call getarg(1, arg(1))
      read(arg(1),*) idoc
      !Sampling size: number of null networks to generate. We recommend a minimum sampling size of 10³ networks.
      call getarg(2, arg(2))
      read (arg(2),*) Nsampl


!****Opening files 
      write (rstr, '(I3.1)') idoc
      generaldoc='general'//trim(adjustl(rstr))//'.txt'
      matrixdoc='matrix'//trim(adjustl(rstr))//'ord.txt'
      mranddoc='matrix'//trim(adjustl(rstr))//'rand.txt'
      resultsdoc='NODF.txt'
      open (10,file=generaldoc,status='old')
      open (11, file=matrixdoc, status='old')
      open(12, file=mranddoc, status='old')

!*****Reading files and allocating variables
      read(10,*) N_R, N_C !Dimension of the bi-adjacency matrix
      N = N_R+N_C !Total number of nodes
      Nprow = (N_R*(N_R-1))/2 !Number of pairs of rows
      Npcol = (N_C*(N_C-1))/2 !Number of pairs of columns
      allocate (M(N_R, N_C), Mrand(N_R, N_C), Prand(N_R, N_C), nr_d(N_R), nc_d(N_C), nodf_vec(Nsampl, 2))
      allocate (nr_drand(N_R), nc_drand(N_C),Npaired_row(Nprow), Npaired_col(Npcol),Mclean(N_R,N_C),Mrand0(N_R,N_C))
      allocate (nr_drandclean(N_R), nc_drandclean(N_C))

!****Reading the real matrix M, the probability matrix Prand, and calculating the degree sequences 
      !Empirical degree sequences of rows in vector nr_d & average null degree sequences of rows in vector nr_drand
      do I=1,N_R
            nr_d(I) = 0
            nr_drand(I) = 0
            read(11,*) (M(I,J), J=1,N_C) !Real bi-adjacency matrix
            read(12,*) (Prand(I,J), J=1,N_C) !Probability bi-adjacency matrix
            do J=1,N_C
                  rbuff = nr_d(I)
                  nr_d(I) = M(I,J) + rbuff
            enddo
      enddo
      !Empirical degree sequences of columns in vector nc_d & average null degree sequences of columns in vector nc_drand
      do J=1,N_C
            nc_d(J) = 0.d0
            nc_drand(J) = 0.d0
            do I=1,N_R
                  rbuff = nc_d(J)
                  nc_d(J) = M(I,J)+ rbuff
            enddo
      enddo


!*******************************************************
!*******CALCULATION OF THE NODF OF THE REAL MATRIX******
!*******************************************************

      !Call of the subroutine that calculates the unnormalised nestedness among columns
      call nested_col(M, nc_d, N_R, N_C, Npcol, Npaired_col)

      !Call of the subroutine that calculates the unnormalised nestedness among rows
      call nested_row(M, nr_d, N_R, N_C, Nprow, Npaired_row)

      !Calculation of the global nestedness
      nest_real = 0.d0
      do I=1, Npcol, 1
            suma = nest_real
            nest_real = suma + Npaired_col(I)/dble(Npcol+Nprow)
      enddo
      do I=1, Nprow, 1
            suma = nest_real
            nest_real = suma + Npaired_row(I)/dble(Npcol+Nprow)
      enddo


!****************************************************************************************
!****CALCULATION OF THE AVERAGE NODF AND ITS STANDARD DEVIATION IN THE NULL ENSEMBLE*****
!****************************************************************************************

!*****Initialising the random number generator
      call dran_ini(123456)

!*****Beggining of the iterations
      do K=1, Nsampl, 1

            !Construction of a new sampled null network
            do I=1,N_R, 1
                  do J=1, N_C, 1
                        u = dran_u()
                        if (u .lt. Prand(I,J))  then
                              Mrand0(I,J) = 1
                        else 
                              Mrand0(I,J) = 0
                        endif
                  enddo
            enddo

            !Ordering the sampled network
            call sormat(N_R,N_C,Mrand0,Mrand)

            !Calculation of the new degree sequences of rows and columns
             do I=1,N_R
                  nr_drand(I) = 0
                  do J=1,N_C
                        nr_drand(I) = Mrand(I,J) + nr_drand(I)
                  enddo
            enddo
            do J=1,N_C
                  nc_drand(J) = 0
                  do I=1,N_R
                        nc_drand(J) = Mrand(I,J) + nc_drand(J)
                  enddo
            enddo

      !*****Cleaning the null network of non-interactive species
            !Counting nodes with degree equal to zero (having no interactions)
            izeros_rows = 0
            do I=N_R,1, -1
                  if (nr_drand(I) .eq. 0) then 
                        izeros_rows = izeros_rows + 1 
                  endif
            enddo
            izeros_col = 0
            do I=N_C,1, -1
                  if (nc_drand(I) .eq. 0) then 
                        izeros_col = izeros_col + 1 
                  endif
            enddo
            !Reallocating the bi-adjacency matrix and removing nodes with no interactions
            deallocate (Mclean)
            allocate (Mclean(N_R-izeros_rows, N_C-izeros_col))
            do I=1, N_R-izeros_rows, 1
                  do J=1, N_C-izeros_col, 1
                        Mclean(I,J) = Mrand(I,J)
                  enddo
            enddo

            !Recalculating the dimension and reallocating vectors for the null network
            N_Rclean = N_R-izeros_rows !Number of rows
            N_Cclean = N_C-izeros_col !Number of columns
            Nprow = (N_Rclean*(N_Rclean-1))/2 !Number of pairs of rows
            Npcol = (N_Cclean*(N_Cclean-1))/2 !Number of pairs of columns
            deallocate(Npaired_col, Npaired_row, nc_drandclean, nr_drandclean)
            allocate (Npaired_row(Nprow), Npaired_col(Npcol), nr_drandclean(N_Rclean), nc_drandclean(N_Cclean))
            do I=1,N_Rclean, 1
                  nr_drandclean(I) =  nr_drand(I)
            enddo
            do I=1,N_Cclean, 1
                  nc_drandclean(I) =  nc_drand(I)
            enddo

     !******Call of the subroutine that calculates the unnormalised nestedness among columns
            call nested_col(Mclean, nc_drandclean, N_Rclean, N_Cclean, Npcol, Npaired_col)

     !******Call of the subroutine that calculates the unnormalised nestedness among rows
            call nested_row(Mclean, nr_drandclean, N_Rclean, N_Cclean, Nprow, Npaired_row)

     !******Calculation of the global nestedness
            nest_rand = 0.d0
            do I=1, Npcol, 1
                  suma = nest_rand
                  nest_rand = suma + Npaired_col(I)/dble((N_R*(N_R-1))/2 + (N_C*(N_C-1))/2)
            enddo
            do I=1, Nprow, 1
                  suma = nest_rand
                  nest_rand = suma + Npaired_row(I)/dble((N_R*(N_R-1))/2 + (N_C*(N_C-1))/2)
            enddo

            nodf_vec(K,1) = nest_rand
            nodf_vec(K,2) = nest_rand*nest_rand

      enddo
      !End of the iterative loop used to sample of the ensemble

      !****Calculation of the average and standard deviation of NODF over the sampling
      suma = 0.d0
      suma2 = 0.d0
      do I=1, Nsampl, 1
            suma = suma + nodf_vec(I,1)
            suma2 = suma2 + nodf_vec(I,2)
      enddo
      nest_rand = suma/dble(Nsampl) 
      sigma = suma2/dble(Nsampl) - nest_rand*nest_rand
      sigma = dsqrt(sigma)


!**************************************************
!******WRITING FINAL RESULTS******
!*************************************************
      open(13, file=resultsdoc, status='unknown')
      write(13,*) "Index of the network=",idoc, "Real NODF=", nest_real, "Average NODF=",nest_rand, "Standard deviation=",sigma

      !Closing files
      close(10)
      close(11)
      close(12)
      close(13)

      end program NODF  
!************************************************************************


!**************RANDOM NUMBERS GENERATOR*****
!******External random number generator 'dranxor.f90', by Raul Toral (see reference [2] above)************
!******Please when using this program ackownledge as well these authors by citing [2]

      INCLUDE 'dranxor.f90'


!****************************************************
!**************SUBROUTINE nested_col****************
!****************************************************
!****This subroutine calculates the unnormalised contribution to NODF of all pairs of columns**** 

      subroutine nested_col(M,MT,N_R,N_C,Npcol,Npaired_col)
      implicit none
      integer :: PO, I, J, K, L, Npcol, N_C, N_R, MT(N_C), M(N_R,N_C)
      real*8 :: Npaired_col(Npcol),  buff, suma 
      K=1
      do I=1, N_C-1, 1
            do J=I+1,N_C,1
                  if (MT(I).gt. MT(J)) then
                        suma = 0.d0
                        do L=1, N_R, 1
                              buff = suma
                              suma = buff + dble(M(L,I))*dble(M(L,J))
                        enddo
                        Npaired_col(K) = 100.d0*dble(suma)/dble(MT(J))
                  elseif (MT(I) .le. MT(J)) then
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

      subroutine nested_row(M,MT,N_R,N_C,Nprow,Npaired_row)
      implicit none
      integer :: PO, I, J, K, L, Nprow, N_R,N_C, MT(N_R), M(N_R,N_C)
      real*8 :: Npaired_row(Nprow), buff, suma 
      K=1
      do I=1, N_R-1, 1
            do J=I+1,N_R,1
                  if (MT(I).gt. MT(J)) then
                        suma = 0
                        do L=1, N_C, 1
                              buff = suma
                              suma = buff + dble(M(I,L))*dble(M(J,L))
                        enddo
                        Npaired_row(K) = 100.d0*dble(suma)/dble(MT(J))
                  elseif (MT(I) .le. MT(J)) then
                        Npaired_row(K)=0
                  endif
                  K=K+1
            enddo 
      enddo
      end


!****************************************************
!**************SUBROUTINE sormat**************
!****************************************************
!****This subroutine orders the bi-adjacency matrix by decreasing degree, both for rows & columns**** 
!****It calls another subroutine called sorvec

      subroutine sormat(nf,nc,mat,matf)
      dimension mat(nf,nc),matf(nf,nc), nvfil(4501),nvcol(4501)
      dimension kufil(4501),kucol(4501), matx(4501,4501) 
      dimension lfil(4501),lcol(4501)

      do 20 i1=1,nf
      nsum=0
      lcol(i1)=0
      do 10 i2=1,nc
      if(mat(i1,i2).ne.0)lcol(i1)=i2
   10 nsum=nsum+mat(i1,i2)
   20 nvfil(i1)=nsum
      call sorvec(nf,nvfil,lcol,kufil)

      do 40 i2=1,nc
      nsum=0
      lfil(i2)=0
      do 30 i1=1,nf
      if(mat(i1,i2).ne.0)lfil(i2)=i1
   30 nsum=nsum+mat(i1,i2)
   40 nvcol(i2)=nsum
      call sorvec(nc,nvcol,lfil,kucol)

      do 60 i1=1,nf
      ix=1+nf-kufil(i1)
      do 50 i2=1,nc
   50 matx(ix,i2)=mat(i1,i2)
   60 continue
      do 80 i2=1,nc
      ix=1+nc-kucol(i2)
      do 70 i1=1,nf
   70 matf(i1,ix)=matx(i1,i2)
   80 continue
      return
      end


!***********SUBROUTINE sorvec**********
! ****************************************
!****Auxiliary subroutine called by sormat, which orders the bi-adjacency matrix by decreasing degree
      subroutine sorvec(nucom,nvin,lfc,kfin)
      dimension nvin(4501),kfin(4501)
      dimension lfc(4501)


      do 100 i1=1,nucom
      ix=1
      do 50 i2=1,nucom
      if(i2.eq.i1)go to 50
      if(nvin(i1).lt.nvin(i2))go to 50
      if(nvin(i1).eq.nvin(i2))go to 30
      ix=ix+1
      go to 50
   30 if(lfc(i1).lt.lfc(i2))go to 50
      if(lfc(i1).eq.lfc(i2))go to 40
      ix=ix+1
      go to 50
   40 if(i2.lt.i1)ix=ix+1

   50 continue
      kfin(i1)=ix
  100 continue
      return
      end
