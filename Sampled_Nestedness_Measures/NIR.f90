      module global
            integer :: N_R, N_C !Global variables: number of rows, number of columns
      end module global

!!*********************************************************************************************
!***********************************PROGRAM: NIR********************************

!     This program calculates, for the nestedness metrics NIR (nesting index based on robustness):
              !i) the real value of nestedness of an empirical network
              !ii) the first two moments of the metric, in a maximum entropy ensemble, computed by numerically sampling the ensemble
!     It takes as input the ordered real bi-adjacency matrix as well as the probability matrix over the null ensemble
!     It returns as output the real measure of NIR, the average expected value of NIR in the ensemble and the standard deviation of NIR in the ensemble

!     Code by Clàudia Payrató Borràs, Laura Hernández and Enrique Burgos. For any questions on the code, you may contact us in the following mail: claudiapb13@gmail.com
!     When using this program, please cite: 
!     [1] "Breaking the Spell of Nestedness: The Entropic Origin of Nestedness in Mutualistic Systems", C. Payrató-Borràs et al. 
!     Phys. Rev. X 9 (2019)
!     [2] "Generation of Gaussian distributed random numbers by using a numerical inversion method, R. Toral, A. Chakrabarti, Computer
      !Physics Communications 74, 327 (1993)

!*********************************************************************************************

      PROGRAM NIR
      use global
      implicit none
      character(len=25) arg(2), generaldoc, matrixdoc, mranddoc, resultsdoc, rstr
      integer, allocatable :: adj(:,:),adj_tmp(:,:),np_dead(:),np_surv(:)
      integer, allocatable :: nr_d(:),nc_d(:),na_dead(:),na_surv(:), Mrand(:,:), Mrand0(:,:) 
      real*8, allocatable :: f_an_killed(:),f_dead_pl(:),f_surv_pl(:), Prand(:,:), nir_vec(:,:)
      real*8, allocatable :: f_pl_killed(:),f_dead_an(:),f_surv_an(:), nc_dran(:), nr_dran(:)
      integer :: vnumber, surv, idoc, I, J, K, ik, iter, Nsamp, iz_r, iz_c, N_Rcn, N_Ccn
      real*8 :: del, den, phi, nest_coeff_pl, nest_coeff_an, A1, A2, nest_an, nest_pl, global_nest
      real*8 :: dran_u, ran, suma, suma2, nest_real, nest_rand, sigma
 

!*******************************************************
!******DEFINITION AND INITIALIZATION of VARIABLES*******
!*******************************************************

!****External arguments (to be entered by terminal)
      !Network index 'idoc', defines the real network upon which the null model is constructed
      call getarg(1, arg(1))
      read(arg(1),*) idoc
      !Sampling size: number of null networks to generate. We recommend a minimum sampling size of 10³ networks.
      call getarg(2, arg(2))
      read (arg(2),*) Nsamp

!*****Opening files 
      write (rstr, '(I3.1)') idoc
      generaldoc='general'//trim(adjustl(rstr))//'.txt'
      matrixdoc='matrix'//trim(adjustl(rstr))//'ord.txt'
      mranddoc='matrix'//trim(adjustl(rstr))//'rand.txt'
      resultsdoc='NIR.txt'
      open (10,file=generaldoc,status='old') 
      open (11, file=matrixdoc, status='old') 
      open(12, file=mranddoc, status='old') 

!*****Reading files and allocating variables
      read(10,*) N_R, N_C !Dimension of the bi-adjacency matrix
      allocate(adj(N_R,N_C),adj_tmp(N_R,N_C),np_dead(N_C), np_surv(N_C), nr_dran(N_R), nc_dran(N_C))
      allocate(na_dead(N_R),na_surv(N_R), nc_d(N_C), nr_d(N_R), Mrand0(N_R, N_C), Mrand(N_R, N_C), Prand(N_R, N_C))
      allocate(f_an_killed(N_C),f_dead_pl(N_C),f_surv_pl(N_C), nir_vec(Nsamp,2))
      allocate(f_pl_killed(N_R),f_dead_an(N_R),f_surv_an(N_R))

!****Reading the real matrix 'adj' and the probability matrix 'Prand', and calculating the density of links 'phi'
      phi = 0.0
      do I=1, N_R, 1 
            read(11,*) (adj(I,J),J=1,N_C)
            read(12,*) (Prand(I,J), J=1,N_C)
            do J=1, N_C, 1
                  phi = phi + float(adj(I,J))
            enddo
      end do 
!     Calculation of the density of links 'phi', and of its inverse (1-phi) 'den'
      phi = phi/float(N_R*N_C)
      den=1-phi


!*******************************************************
!*******CALCULATION OF THE NIR OF THE REAL MATRIX******
!*******************************************************

!******** NIR for ROWS *******

!*****Node removal by decreasing degree order (DDR strategy)

      !Initialising variables
      adj_tmp=adj !temporal matrix
      np_dead=0
      np_surv=0
      f_an_killed=0.0
      f_dead_pl=0.0
      f_surv_pl=0.0

      !Calculation of the fraction of surviving row nodes when deleting column nodes with DDR strategy
      do i=1,N_C 
            call an_killer(N_R,N_C,adj_tmp,i,vnumber,surv)
            np_dead(i)=vnumber
            np_surv(i)=surv
            f_an_killed(i)=dfloat(i)/dfloat(N_C) 
            f_dead_pl(i)=dfloat(np_dead(i))/dfloat(N_R)
            f_surv_pl(i)=dfloat(np_surv(i))/dfloat(N_R)
      end do

      !Area calculation under this ATC (attack tolerance curve), saved as 'A2'
      del=1./dfloat(N_C)
      call area(N_C,f_surv_pl,del,A2)

!****Node removal by increasing degree order (IDR strategy)

      !Initialising variables
      adj_tmp=adj
      np_dead=0
      np_surv=0
      f_an_killed=0.0
      f_dead_pl=0.0
      f_surv_pl=0.0

      !Calculation of the fraction of surviving row nodes when deleting column nodes with IDR strategy
      do i=1,N_C
            ik=N_C-i+1 
            call an_killer(N_R,N_C,adj_tmp,ik,vnumber,surv)
            np_dead(i)=vnumber
            np_surv(i)=surv
            f_an_killed(i)=dfloat(i)/dfloat(N_C) 
            f_dead_pl(i)=dfloat(np_dead(i))/dfloat(N_R)
            f_surv_pl(i)=dfloat(np_surv(i))/dfloat(N_R)
      end do

      !Calculation of the area under this ATC (attack tolerance curve), saved as 'A1'
      del=1./dfloat(N_C)
      call area(N_C,f_surv_pl,del,A1)

!*****Total NIR for rows
      nest_coeff_pl=A1-A2
      nest_coeff_pl=nest_coeff_pl/den


!******* NIR for COLUMNS *******

!*****Node removal by decreasing degree order (DDR strategy)

      !Initialising variables
      adj_tmp=adj
      na_dead=0
      na_surv=0
      f_pl_killed=0.0
      f_dead_an=0.0
      f_surv_an=0.0

      !Calculation of the fraction of surviving column nodes when deleting row nodes with DDR
      do i=1,N_R 
             call pl_killer(N_R,N_C,adj_tmp,i,vnumber,surv)     
             na_dead(i)=vnumber
             na_surv(i)=surv
             f_pl_killed(i)=dfloat(i)/dfloat(N_R) 
             f_dead_an(i)=dfloat(na_dead(i))/dfloat(N_C)
             f_surv_an(i)=dfloat(na_surv(i))/dfloat(N_C)
      end do

      !Calculation of the area under this ATC (attack tolerance curve), saved as 'A2'
      del=1./dfloat(N_R)
      call area(N_R,f_surv_an,del,A2)

!****Node removal by increasing degree order (IDR strategy)

      !Initialising variables
      adj_tmp=adj
      na_dead=0
      na_surv=0
      f_pl_killed=0.0
      f_dead_an=0.0
      f_surv_an=0.0

      !Calculation of the fraction of surviving column nodes when deleting row nodes with IDR
      do i=1,N_R
               ik=N_R-i+1 
               call pl_killer(N_R,N_C,adj_tmp,ik,vnumber,surv) 
               na_dead(i)=vnumber
               na_surv(i)=surv
               f_pl_killed(i)=dfloat(i)/dfloat(N_R) 
               f_dead_an(i)=dfloat(na_dead(i))/dfloat(N_C)
               f_surv_an(i)=dfloat(na_surv(i))/dfloat(N_C)
      end do

     !Calculation under this ATC (attack tolerance curve), saved as 'A1'
      del=1./dfloat(N_R)
      call area(N_R,f_surv_an,del,A1)
 
!****Total NIR for columns
      nest_coeff_an=A1-A2
      nest_coeff_an=nest_coeff_an/den


!***** GLOBAL NIR for the REAL network **** 
      !This takes into account both rows' and columns' contribution
      nest_real = (float(N_C)*nest_coeff_an+float(N_R)*nest_coeff_pl)/dfloat(N_R+N_C)


!****************************************************************************************
!****CALCULATION OF THE AVERAGE NIR AND ITS STANDARD DEVIATION IN THE NULL ENSEMBLE*****
!****************************************************************************************

!*****Initialising the random number generator
      call dran_ini(123456)

!*****Beggining of the iterations to sample the ensemble
      do iter=1, Nsamp, 1 

            !Construction of a new sampled null network
            do I=1, N_R, 1
                  do J=1, N_C, 1
                        ran = dran_u()
                        if (ran .lt. Prand(I,J))  then
                              Mrand0(I,J) = 1
                        else 
                              Mrand0(I,J) = 0
                        endif
                  enddo
            enddo

            !Ordering the sampled network
            call sormat(N_R,N_C,Mrand0,Mrand0)

            !Calculation of the new degree sequences of rows and columns
             do I=1,N_R
                  nr_dran(I) = 0
                  do J=1,N_C
                        nr_dran(I) = Mrand0(I,J) + nr_dran(I)
                  enddo
            enddo
            do J=1,N_C
                  nc_dran(J) = 0
                  do I=1,N_R
                        nc_dran(J) = Mrand0(I,J) + nc_dran(J)
                  enddo
            enddo

      !*****Cleaning the null network of non-interactive species

            !Counting nodes with degree equal to zero (having no interactions)
            iz_r = 0 !rows
            do I=N_R,1, -1
                  if (nr_dran(I) .eq. 0) then 
                        iz_r = iz_r + 1 
                  endif
            enddo
            iz_c = 0 !columns
            do I=N_C,1, -1
                  if (nc_dran(I) .eq. 0) then 
                        iz_c = iz_c + 1 
                  endif
            enddo

            !Reallocating the bi-adjacency matrix and removing nodes with no interactions
            deallocate (Mrand, adj, adj_tmp)
            allocate (Mrand(N_R-iz_r, N_C-iz_c), adj(N_R-iz_r, N_C-iz_c), adj_tmp(N_R-iz_r, N_C-iz_c))
            do I=1, N_R-iz_r, 1
                  do J=1, N_C-iz_c, 1
                        Mrand(I,J) = Mrand0(I,J)
                  enddo
            enddo
            adj = Mrand 

            !Recalculating the dimension and reallocating vectors for the null network
            N_Rcn = N_R-iz_r !New number of rows
            N_Ccn = N_C-iz_c !New number of columns

            !Recalculating the density of links of the network, 'phi'
            do I=1, N_Rcn, 1
                  do J=1, N_Ccn, 1
                        phi = phi + float(Mrand(I,J))
                  enddo
            end do 
            phi = phi/float(N_Rcn*N_Ccn)
            den=1-phi


      !********NIR for ROWS*******

      !*****Node removal by decreasing degree order (DDR strategy)

            !Initialising variables
            adj_tmp=adj 
            np_dead=0
            np_surv=0
            f_an_killed=0.0
            f_dead_pl=0.0
            f_surv_pl=0.0

            !Calculation of the fraction of surviving row nodes when deleting column nodes with DDR
            do i=1,N_Ccn 
                  call an_killer(N_Rcn,N_Ccn,adj_tmp,i,vnumber,surv)
                  np_dead(i)=vnumber
                  np_surv(i)=surv
                  f_an_killed(i)=dfloat(i)/dfloat(N_Ccn) 
                  f_dead_pl(i)=dfloat(np_dead(i))/dfloat(N_Rcn)
                  f_surv_pl(i)=dfloat(np_surv(i))/dfloat(N_Rcn)
            end do

            !Area calculation under this ATC (attack tolerance curve), saved as 'A2'
            del=1./dfloat(N_Ccn)
            call area(N_Ccn,f_surv_pl,del,A2)


      !*****Node removal by increasing degree order (IDR strategy)

            !Initialising variables
            adj_tmp=adj
            np_dead=0
            np_surv=0
            f_an_killed=0.0
            f_dead_pl=0.0
            f_surv_pl=0.0

            !Calculation of the fraction of surviving row nodes when deleting column nodes with IDR
            do i=1,N_Ccn
                  ik=N_Ccn-i+1 
                  call an_killer(N_Rcn,N_Ccn,adj_tmp,ik,vnumber,surv)
                  np_dead(i)=vnumber
                  np_surv(i)=surv
                  f_an_killed(i)=dfloat(i)/dfloat(N_Ccn) 
                  f_dead_pl(i)=dfloat(np_dead(i))/dfloat(N_Rcn)
                  f_surv_pl(i)=dfloat(np_surv(i))/dfloat(N_Rcn)
            end do

            !Area calculation under this ATC (attack tolerance curve), saved as 'A1'
            del=1./dfloat(N_Ccn)
            call area(N_Ccn,f_surv_pl,del,A1)
      
      !*****Total NIR for rows
            nest_coeff_pl=A1-A2      
            nest_coeff_pl=nest_coeff_pl/den


      !*******NIR for COLUMNS*******

      !*****Node removal by decreasing degree order (DDR strategy)

            !Initializing variables
            adj_tmp=adj
            na_dead=0
            na_surv=0
            f_pl_killed=0.0
            f_dead_an=0.0
            f_surv_an=0.0

            !Calculation of the fraction of surviving column nodes when deleting row nodes with DDR
            do i=1,N_Rcn
                   call pl_killer(N_Rcn,N_Ccn,adj_tmp,i,vnumber,surv)     
                   na_dead(i)=vnumber
                   na_surv(i)=surv
                   f_pl_killed(i)=dfloat(i)/dfloat(N_Rcn) 
                   f_dead_an(i)=dfloat(na_dead(i))/dfloat(N_Ccn)
                   f_surv_an(i)=dfloat(na_surv(i))/dfloat(N_Ccn)
            end do

            !Calculation of the area under this ATC (attack tolerance curve), saved as 'A2'
            del=1./dfloat(N_Rcn)
            call area(N_Rcn,f_surv_an,del,A2)

      !*****Node removal by increasing degree order (IDR strategy)

            !Initializing variables
            adj_tmp=adj
            na_dead=0
            na_surv=0
            f_pl_killed=0.0
            f_dead_an=0.0
            f_surv_an=0.0

            !Calculation of the fraction of surviving column nodes when deleting row nodes with IDR
            do i=1,N_Rcn
                     ik=N_Rcn-i+1 
                     call pl_killer(N_Rcn,N_Ccn,adj_tmp,ik,vnumber,surv) 
                     na_dead(i)=vnumber
                     na_surv(i)=surv
                     f_pl_killed(i)=dfloat(i)/dfloat(N_Rcn) 
                     f_dead_an(i)=dfloat(na_dead(i))/dfloat(N_Ccn)
                     f_surv_an(i)=dfloat(na_surv(i))/dfloat(N_Ccn)
            end do

            !Calculation of the area under this ATC (attack tolerance curve), saved as 'A1'
            del=1./dfloat(N_Rcn)
            call area(N_Rcn,f_surv_an,del,A1)
      
      !*****Total NIR for columns
            nest_coeff_an=A1-A2
            nest_coeff_an=nest_coeff_an/den


      !*****GLOBAL NIR for the sampled null network (taking into account both rows and columns)
            global_nest=(float(N_Ccn)*nest_coeff_an+float(N_Rcn)*nest_coeff_pl)/dfloat(N_Rcn+N_Ccn)
            nir_vec(iter,1) = global_nest 
            nir_vec(iter,2) = global_nest*global_nest

      enddo
      !End of the iterative loop used to sample of the ensemble


      !****Calculation of the average and standard deviation of NIR over the sampling
      suma = 0.d0
      suma2 = 0.d0
      do I=1, Nsamp, 1
            suma = suma + nir_vec(I,1)
            suma2 = suma2 + nir_vec(I,2)
      enddo
      nest_rand = suma/dble(Nsamp) !average
      sigma = suma2/dble(Nsamp) - nest_rand*nest_rand 
      sigma = dsqrt(sigma) !standard deviation


!**************************************************
!******WRITING FINAL RESULTS******
!*************************************************
      open(13, file=resultsdoc, status='unknown')
      write(13,*) "Index of the network=",idoc, "Real NIR=", nest_real, "Average NIR=",nest_rand, "Standard deviation=",sigma

      !Closing files
      close(10)
      close(11)
      close(12)
      close(13)

      end program NIR
!************************************************************************


!**************RANDOM NUMBERS GENERATOR*****
!******External random number generator 'dranxor.f90', by Raul Toral (see reference [2] above)************
!******Please when using this program ackownledge as well these authors by citing [2]

      INCLUDE 'dranxor.f90'

!********************************************************************** 
        subroutine an_killer (N_R,N_C,adj_tmp,ik,vnumber,surv)
!**********************************************************************
! Kills nodes from columns and modifies an evolutive adjacency matrix erasing the corresponding nodes from columns and the derived dead nodes from rows.
 
        implicit real*8 (a-h,o-z)
        integer  vnumber,ik,N_R,N_C, surv
        integer :: adj_tmp(N_R, N_C)

        vnumber=0
        surv=0
        do i=1,N_R
           adj_tmp(i,ik)=0
        end do
        do i=1,N_R
           ktmp=0
           do j=1,N_C
              ktmp=ktmp+adj_tmp(i,j)
           end do
              if (ktmp==0)then
               vnumber=vnumber+1
              else
                surv=surv+1
            end if
         end do
           return
           end
       
!********************************************************************** 
        subroutine pl_killer (N_R,N_C,adj_tmp,ik,vnumber,surv)
!**********************************************************************
! Kills nodes from rows and modifies an evolutive adjacency matrix erasing the corresponding rows and the derived dead nodes from columns.

        implicit real*8 (a-h,o-z)
        integer  vnumber,ik,nr,N_R,N_C, surv
        integer :: adj_tmp(N_R,N_C)

        vnumber=0
        surv=0
        do i=1,N_C
           adj_tmp(ik,i)=0
        end do
        do j=1,N_C
           ktmp=0
           do i=1,N_R
              ktmp=ktmp+adj_tmp(i,j)
           end do
              if (ktmp==0)then
               vnumber=vnumber+1
              else
                surv=surv+1
	      end if
         end do
           return
           end

!************************************************                      
         subroutine area(nspecies,fsurv,del,A)
!************************************************
         integer nspecies
         double precision, dimension (nspecies):: fsurv
         double precision A,del
         A=0
          do i=1,nspecies
            A=A+fsurv(i)
         end do

         A=A*del
         return
         end    
!**************************************************


!***********SUBROUTINE sormat**********
!*******reorders the matrix ************
      subroutine sormat(nf,nc,mat,matf)
      dimension mat(nf,nc),matf(nf,nc), nvfil(4501),nvcol(4501)
      dimension kufil(4501),kucol(4501), matx(4501,4501) 
      dimension lfil(4501),lcol(4501)

!      kfil=0
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

!  las filas y columnas han sido numeradas en orden creciente del numero de
!  interacciones. En las matrices triang.-inf.-ezq. tal orden es decreciente

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
!*************************************

!***********SUBROUTINE sorvec**********
! ****************************************
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
!***************************************






