!!!!!!!!!!!!!!!!!!!!!!!
!!! Modification: TONI (LHref)
!!! This subroutine selects a random LH bead that can be moved in the elastic model
!!!
!!! Input
!!!  - nCtNt: number of N-term and C-term beads in a LH
!!!  - nLHmove: total number of LH beads that can be selected (all nucleosomes)
!!!  - nbN: number of N-term beads
!!!  - nbG: number of globular-head beads
!!!  - nbH: number of beads in one LH
!!!  - LHb: selected beads within the range h_n
!!!!!!!!!!!!!!!!!!!!!!!
      subroutine LHbead(nCtNt,nLHmove,nbN,nbG,nbH,LHb,np,myid,ierr)



      use modglob
      use mpi
      implicit NONE
  
  
  !!! Variables

      integer nCtNt,nbN,nbG,nbH
      integer nLHmove,LHb
      integer i,j,k
      integer ic,ib,ih
  	  integer flbound
  	  integer myid,ierr,np
      double precision random1,ranu


      if (myid.eq.0) then
       random1 = ranu()
      endif
  
	  call MPI_BCAST(random1, 1, MPI_DOUBLE_PRECISION, 
     +	    0,MPI_COMM_WORLD,ierr)
     
C      write (*,*) 'INSIDE LHBEAD'

      !!!!!!!!!!!!!!!!!!!!!
      !!! Modification: TONI (LH concentration)
      !!! Peak a bead from a LH that is bound to a nucleosome
      flbound = 0
      do while(flbound.eq.0)
      !!!!!!!!!!!!!!!!!!!!!

      !!! Select a random bead in the range nLHmove
      i = (nLHmove*random1) + 1
      if(i.gt.nLHmove) i = nLHmove
  
      !!! Find the bead in the h_n range
      ic = (i-1)/nCtNt + 1 !!! associated core
      ib = mod(i,nCtNt) !!! order within the LH Ct-Nt
      if(ib.eq.0) ib = nCtNt
      if(ib.LE.nbN) then 
        ih = (ic-1)*nbH + ib !!! bead in the h_n range (N-term)
      else
        ih = (ic-1)*nbH + ib + nbG !!! bead in the h_n range (C-term)
      endif
   !write(*,*) 'i,ic,ib',i,ic,ib,ih

   !!! Bead
      LHb = ih

  !!!!!!!!!!!!!!!!!!!!!
  !!! Modification: TONI (LH concentration)
      if(LHboundbds(LHb)) then
        flbound = 1
      endif
      enddo !!! while(flbound.eq.0)
  !!!!!!!!!!!!!!!!!!!!!


C      write(*,*) 'FINISHED LHBEAD'
      return
      end subroutine LHbead
