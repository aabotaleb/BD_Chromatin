ccccccccccccccc
c	  Mesoscale Metropolis MC Modeling of Chromatin Fibers as 
C	  remodeled by G.BASCOM ~10/14	
c     This code can simulate single fibers with different LH concentrations.
c	  Makes use of DiSCO for nucleosomes, beads for DNA, beads for Linker Histones
c	  stretch, bend and torsion done in the usual manner
c	  Excluded volume by 12-6 Lennerd Jones
c	  Long Range electrostatics by Debye-Huckel
c	  Tails regrown by Rosenbluth and Rosenbluth
cccccccccccccccc

      program main
      use modglob
      use mpi
      implicit NONE
	  integer n_c,nlb,nc3,n,n3,t_n,t_n3,Nq,Nq3,h_n,h_n3,np,ierr,myid
CCCCCCC^_^_^_^_^_^^_^_^_^_^_^^_^_^_^_^_^^_^_^_^_^_^^_^_^_^_^_^CCCCCCCCC
CCCCCCC^_^_^_^_^_^^_^_^_^_^_^^_^_^_^_^_^^_^_^_^_^_^^_^_^_^_^_^CCCCCCCCC        
CCCCCCC change vector nbn length to 750 rather then 200 as n_c in dim.in read equals 750
        integer A,nbn(750),i,nlbtot 
	  character(LEN=10) dim_file
      real*4 tini,tend
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCC Debugging Purposes  Abdelrhman Modifications 8 Sept 2015         CCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Commenting the following lines:                       
C       A = IARGC()
C	  IF (A < 4) THEN
C	  write(*,*) 'ERROR: WRONG NUM OF ARGS'
C	  write(*,*) 'USE./exe dim_file para_file LH_in LH_equil > out'
C	  GO TO 9990
C	  ENDIF
c and  call GETARG(1, dim_file)  replaced by   dim_file='dim.in'
c and  call GETARG(2, inputfile) replaced by   inputfile='input.run' 
c and  call GETARG(3,LH_input)   replaced by   LH_input='LH.in'
c and  call GETARG(4,LH_equil)   replaced by   LH_equil='LH_equil.in'

       
C 	Open dim_file to get the dimensions with format:
C 	number_of_cores, number_of_beads
C
C 	beads_in_NTERM,beads_in_Globular_Head,beads_in_CTERM

      call cpu_time(tini)
      
      
C 	Initialize MPI
       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)

      if (myid.eq.0) then
      print*,'**************************************************
     +********************'
      print*, "*   Mechnical Forces of Chromatin Mesoscale Model  
     +Calculations   *" 
      print*,'*                  Analytical Version Only with TestGH    
     +           *'
      print*,'**************************************************
     +********************'
      endif
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCC Debugging Purposes  Abdelrhman Modifications 8 Sept 2015         CCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C            call GETARG(1, dim_file) is replaced by dim_file='dim.in'
C     call GETARG(1, dim_file)
      dim_file='t1/dim.in'
       
      open(11,file=dim_file)
       write(*,'(A,i2)'),' DIMENSION FILE BEING READ BY.',myid
      read(11,*)
      read(11,*)n_c
      read(11,*)
      read(11,*)nbNt,nbGh,nbCt
      read(11,*)
      
C	  **pulled from rosana's vNRL code, read in linker bead lengths'      
      
      nlbtot=0
      do i=1,n_c
         read(11,*)nbn(i)
         nlbtot=nlbtot+nbn(i)
      enddo
      
      
      
      nbLH = nbNt + nbGh + nbCt
      close(11)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      


c	Calculate major vars for allocation:

       nc3=n_c*3 		!number of coresx3 
       n=n_c+nlbtot 	!total number of dna + core beads 
       n3=3*n 			!total number of coors
       t_n=n_c*50 		!total number of histone tail beads
       t_n3=3*t_n         !total number of coors for tail beads -Abotaleb- 
       Nq = 300 		!total number of DiSCO charges
       Nq3 = 3*Nq		!total number of DiSCO charge coords
       h_n = n_c*nbLH	+1!total number of LH beads
	   h_n3=3*h_n         !total number of coordinates for LH beads -Abotaleb-


       
C       write (*,*) 'total num of dna beads', n
C       write (*,*) 'total num of disco charges', Nq
C       write (*,*) 'total num of LH beads', h_n
C       write (*,*) 'total num of tail beads', t_n



C 	Run Main Program


		
       call init(nlb,n_c,nc3,n,n3,t_n,t_n3,Nq,Nq3,
     + h_n,h_n3,ierr,myid,np,nbn)
		
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)		
	   if (myid.eq.0) then
       WRITE(*,*)"******************************"
       WRITE(*,*)"* PROGRAM COMPLETED NORMALLY *"
       WRITE(*,*)"******************************"

       call cpu_time(tend)
       WRITE(*,*) 'CPU TIME ELAPSED:',(tend-tini)/60,'MINS'
C     +  									  (tend-tini),'SECONDS'
       endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
 9990  end
     




     	
