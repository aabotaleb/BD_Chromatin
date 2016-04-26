!!!!!!!!!!!!!!!!!!!!!!!
!!! Modification: TONI (LH concentration)
!!! This subroutine distributes the LHs in nucleosmes 
!!! for given concentration or specific input distribution
!!!
!!! Inputs
!!!  - myid: id associated to the trajectory
!!!  - n_c: number of nucleosomes
!!!  - h_n: number of LH beads
!!!  - seed: seed for the random number generator
!!!  - LHconc: concentration of LHs
!!!  - modeLHc: distribution mode (random or input)
!!!  - LHnum: number of LHs (necessary for modeLHc=2)
!!!
!!! Outputs
!!!  - LHbound : List of LHs (bound=1,unbound=0)
!!!  - LHboundbds : List of LH beads (bound=1,unbound=0)
!!!!!!!!!!!!!!!!!!!!!!!
subroutine LHboundsub(myid,np,n_c,h_n,seed,LHconc,modeLHc,LHnum,ierr,VERBOSE)

  !!!   We make use of global parameters/variables
  ! LHbound
  ! LHboundbds
  use modglob
  use mpi


  !!! Variables
  implicit NONE
  integer jmin,jmax,kmin,kmax
  double precision skip
  double precision LHconc  !!! Concentration of LHs/nucleosome
  integer modeLHc          !!! Mode for the LH concentration algorithm
  integer LHcount,LHnum,np
  integer n_c,h_n
  integer i,j,k
  integer seed
  character ifile*50, ofile*50
  integer myid
  double precision PLHn
  integer flbound
  integer ierr,z
  logical VERBOSE


  !!! Output file name
  ofile = 'LHbound.0.out'

  !!! Input file name
  ifile = 'LHbound.0.in'



  !!! Allocate memory
  allocate(LHbound(n_c))    !!! Bound vector per LH (bound=1, unbound=0)
  allocate(LHboundbds(h_n)) !!! Bound vector per LH bead

! if (myid.eq.0).and.(VERBOSE) then
!   do i = 1,h_n
!    print*,LHbound_randarray(i)
!   enddo
! endif
! ...................................................START MASTER LOOP
!  if (myid.eq.0) then


  !!! Starting random number generator
  !!! Notice that this generator is different than in the main program
  !!! In this way, we don't alter the sequence of random numbers of ranu()
!  call random_seed(put=seed)


  !!! Assign LHs to nucleosomes
  if(modeLHc.eq.0) then

     do i=1,n_c
        if(LHbound_randarray(i).le.LHconc) then
           LHbound(i) = 1
        else
           LHbound(i) = 0
        endif
     enddo
       	 
        

  elseif(modeLHc.eq.1) then
     !!! Manual distribution (input file)
     !open(unit=101,file="LHbound.in")
     open(unit=101, file=ifile)
     print*,'READING LHbound.in on proc',myid
     read(101,*)
     read(101,*)
     do i=1,n_c
        read(101,*) LHbound(i)
     enddo

  else
     !!! Random distribution of a specific number of LHs
     LHbound(:) = 0  
     do i=1,LHnum
     !!! Randomly pick a nucleosome without LH (from 1 to n_c)
        flbound = 0
        z=0
        do while(flbound.eq.0)
!           call random_number(Pdice)
           j = nint((n_c - 1)*LHbound_randarray(z) + 1) !!! random core between 1 and n_c
           z=z+1
           if(LHbound(j).eq.0) then
              LHbound(j) = 1           
              flbound = 1
           endif        
        enddo
     enddo
  endif
  
      if ((myid.eq.0).and.(VERBOSE)) then
        print*,'The Cores Marked (1) will have LHs bound:'
      		do i=1,n_c
       			print*,'CORE',i,':',LHbound(i) 
      		end do
	   endif


  !!! Distribution of LH beads | counting LHs
  LHcount = 0
  do i=1,n_c
     jmin = (i-1)*nbLH + 1
     jmax = jmin + (nbLH - 1)
     LHboundbds(jmin:jmax) = LHbound(i)

     if(LHbound(i)) then
        LHcount = LHcount + 1             
     endif
  enddo




  !!! Output file
  if (myid.eq.0) then
   open(unit=101, file=ofile)
   write(101,*) '# [LH] (distr) = ',LHconc,'   mode = ',modeLHc
   write(101,*) '# [LH] (fiber) = ',LHcount/(1.*n_c),'   LHs  = ',LHcount
   do i=1,n_c
      write(101,*) LHbound(i)
   enddo
   close(101)
   
  endif

!.......................................................END MASTER LOOP

!  end if 
!!  Broadcast all the variables to other processes
!  LHcount
!  LHbound
!  LHboundbds

!	  call MPI_BCAST(LHcount, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!	  call MPI_BCAST(LHbound, n_c, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!	  call MPI_BCAST(LHboundbds, h_n, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	  
!	  write (*,*) 'LHcount, LHbound, LHboundbds on proc:',myid,LHcount,LHbound,LHboundbds
!
!  do i=1,n_c
!     write(*,*) LHbound(i)
!  enddo



  return
end subroutine LHboundsub
