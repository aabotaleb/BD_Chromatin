cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine global(n_c,n,n3,type,r,a,b,c,r_n,a_n,b_n,c_n,ro,t_n, 
     +     t_X_n,t_Y_n,t_Z_n,t_X,t_Y,t_Z,h_n,h_X_n,h_Y_n,h_Z_n,h_X,h_Y,
     +     h_Z,delta,first,last,pi,h_X0,h_Y0,h_Z0,
     +	   randarray,np,myid,ierr,withlink)


      use modglob
	  use mpi
      implicit NONE

      integer 		n_c,n,n3,type(n),n_b
      integer 		i,i1,i2,i3,j,m1,m2,m3
      integer 		k,ki,k1,k2,k3,l,init
      integer 		first,last, t_n,h_n
      integer         z,np,myid,ierr
      
      double precision	ranu,magor,magcom,maga,magb,magc
      double precision  x1,y1,z1,dist1,dist2,norm,invnorm
      double precision  costh,sinth,phi,theta,delta,mag      
      double precision 	r(n3),a(n3),b(n3),c(n3),ro,pi
      double precision 	r_n(n3),a_n(n3),b_n(n3),c_n(n3)
      double precision 	rp(3),rc(3),rm(3),rmp(3)
      double precision 	par(3),perp(3),pp(3),ss(3)
      double precision 	t_X_n(t_n), t_Y_n(t_n), t_Z_n(t_n)
      double precision 	t_X(t_n), t_Y(t_n), t_Z(t_n)
      double precision 	h_X_n(h_n), h_Y_n(h_n), h_Z_n(h_n)
      double precision 	h_X(h_n), h_Y(h_n), h_Z(h_n)
      double precision  comp(3),h_X0(3),h_Y0(3),h_Z0(3)
      double precision  randarray(17)
      
	double precision mAc(3,3),mAcn(3,3),mAcninv(3,3)
	double precision vcomp(3),txyz(3),rxyz(3),drtxyz(3)
	double precision vrot(3),rprn(3)
	double precision mB(3,3),tnxyz(3)
	double precision ac(3),bc(3),cc(3)
	double precision ac_n(3),bc_n(3),cc_n(3)
      double precision rnxyz(3),proj
      
      logical withlink
      
      interface
         function finvM3x3(mA)
           double precision,dimension(3,3)::mA,finvM3x3
         end function finvM3x3
      end interface
	
      integer index,indexj(n)

	double precision xrel,yrel,zrel


c get the random numbers on the main proc
c broadcast them to the other processors

C      if (myid.eq.0) then
C	    do z = 1, 5 
C         random(z) = ranu()
C	    enddo
C	  endif
C	  
C	  call MPI_BCAST(random, 5, MPI_DOUBLE_PRECISION, 
C     +	    0,MPI_COMM_WORLD,ierr)


c     choose the system component about which to rotate the smallest end of
c     the oligonucleosome about 

c      if(init.eq.0).and.(type(n_b).eq.0)then
c      magor=dsqrt((h_X(ki)-r(m1))**2+(h_Y(ki)-r(m2))**2+
c      +        (h_Z(ki)-r(m3))**2)
c      init=1
c      endif

      index=0
      do i=1,n
         if(type(i).ne.0) index=index+1
         indexj(i)=index
      enddo



      n_b = n*randarray(4) + 1


      if(n_b.gt.n) n_b = n
      if (n_b .ge. n/2) then
         first = n_b
         last = n
      else
         first = 1
         last = n_b - 1
         if (type(n_b) .ne. 0) last = n_b
      endif
     

      costh = 2.0d0*(0.5d0-randarray(5))
      
      if(randarray(6).gt.0.5d0)then
      sinth = dsqrt(1.0d0-costh*costh)
      else
      sinth = -1.d0*dsqrt(1.0d0-costh*costh)
      endif

      phi = randarray(7)*2.d0*pi
      rp(1) = sinth*dcos(phi)
      rp(2) = sinth*dsin(phi)
      rp(3) = costh

ccccccccccc
ccc Modification: DURBA and TONI
c     rp is the random axis where the rotation will be performed.
c     It doesn't require any amplitud of rotation (all space)
ccccccccccccccc'

      theta = delta*(randarray(8) - 0.5d0)
       

      sinth = dsin(theta)
      costh = dcos(theta)
      rc(1) = r(3*(n_b-1)+1)
      rc(2) = r(3*(n_b-1)+2)
      rc(3) = r(3*(n_b-1)+3)

! if core:
      if (type(n_b) .ne. 0) then
          rc(1) = rc(1) - ro*b(3*(n_b-1)+1)
          rc(2) = rc(2) - ro*b(3*(n_b-1)+2)
          rc(3) = rc(3) - ro*b(3*(n_b-1)+3)
      end if

      do k = first, last

          k1 = 3*(k-1)+1
          k2 = 3*(k-1)+2
          k3 = 3*(k-1)+3
!	  
c (1) Translation of coordinates
          rm(1) = r(k1) - rc(1)
          rm(2) = r(k2) - rc(2)
          rm(3) = r(k3) - rc(3)
          mag = rp(1)*rm(1) + rp(2)*rm(2) + rp(3)*rm(3)
          rmp(1) = rm(1) - mag*rp(1)
          rmp(2) = rm(2) - mag*rp(2)
          rmp(3) = rm(3) - mag*rp(3)
          ss(1) = rp(2)*rmp(3) - rp(3)*rmp(2)
          ss(2) = rp(3)*rmp(1) - rp(1)*rmp(3)
          ss(3) = rp(1)*rmp(2) - rp(2)*rmp(1)
          r_n(k1) = rc(1) + costh*rmp(1) + sinth*ss(1) + mag*rp(1)
          r_n(k2) = rc(2) + costh*rmp(2) + sinth*ss(2) + mag*rp(2)
          r_n(k3) = rc(3) + costh*rmp(3) + sinth*ss(3) + mag*rp(3)
c
c (2) Movement of local coordinate systems, {u,f,v}.
c     (a) moving a-vector
          mag = a(k1)*rp(1) + a(k2)*rp(2) + a(k3)*rp(3)
          par(1) = mag*rp(1)
          par(2) = mag*rp(2)
          par(3) = mag*rp(3)
          perp(1) = a(k1) - par(1)
          perp(2) = a(k2) - par(2)
          perp(3) = a(k3) - par(3)
          pp(1) = ( perp(2)*rp(3) - perp(3)*rp(2) )
          pp(2) = ( perp(3)*rp(1) - perp(1)*rp(3) )
          pp(3) = ( perp(1)*rp(2) - perp(2)*rp(1) )
          a_n(k1) = par(1) + costh*perp(1) - sinth*pp(1)
          a_n(k2) = par(2) + costh*perp(2) - sinth*pp(2)
          a_n(k3) = par(3) + costh*perp(3) - sinth*pp(3)
! renormalize
          norm=dsqrt(a_n(k1)**2+a_n(k2)**2+a_n(k3)**2)
          invnorm=1.d0/norm
          a_n(k1)=a_n(k1)*invnorm
          a_n(k2)=a_n(k2)*invnorm
          a_n(k3)=a_n(k3)*invnorm
    
c
c     (b) moving b-vector
          mag = b(k1)*rp(1) + b(k2)*rp(2) + b(k3)*rp(3)
          par(1) = mag*rp(1)
          par(2) = mag*rp(2)
          par(3) = mag*rp(3)
          perp(1) = b(k1) - par(1)
          perp(2) = b(k2) - par(2)
          perp(3) = b(k3) - par(3)
          pp(1) = ( perp(2)*rp(3) - perp(3)*rp(2) )
          pp(2) = ( perp(3)*rp(1) - perp(1)*rp(3) )
          pp(3) = ( perp(1)*rp(2) - perp(2)*rp(1) )
          b_n(k1) = par(1) + costh*perp(1) - sinth*pp(1)
          b_n(k2) = par(2) + costh*perp(2) - sinth*pp(2)
          b_n(k3) = par(3) + costh*perp(3) - sinth*pp(3)

ccccccccccccc
ccc Modification: DURBA and TONI
c     We correct orthogonality of b_n based on a_n (Gram-Schmidt)
          proj = a_n(k1)*b_n(k1)+a_n(k2)*b_n(k2)+a_n(k3)*b_n(k3)
          b_n(k1) = b_n(k1) - proj*a_n(k1)
          b_n(k2) = b_n(k2) - proj*a_n(k2)
          b_n(k3) = b_n(k3) - proj*a_n(k3)
ccccccccccccc
! renormalize
          norm=dsqrt(b_n(k1)**2+b_n(k2)**2+b_n(k3)**2)
          invnorm=1.d0/norm
          b_n(k1)=b_n(k1)*invnorm
          b_n(k2)=b_n(k2)*invnorm
          b_n(k3)=b_n(k3)*invnorm

c     (c) new c-vector
          c_n(k1) = a_n(k2)*b_n(k3) - a_n(k3)*b_n(k2)
          c_n(k2) = a_n(k3)*b_n(k1) - a_n(k1)*b_n(k3)
          c_n(k3) = a_n(k1)*b_n(k2) - a_n(k2)*b_n(k1)
	  
          norm=dsqrt(c_n(k1)**2+c_n(k2)**2+c_n(k3)**2)
          invnorm=1.d0/norm
          c_n(k1)=c_n(k1)*invnorm
          c_n(k2)=c_n(k2)*invnorm
          c_n(k3)=c_n(k3)*invnorm

      enddo
      

c update tail and linker histone positions of the relevant nucleosomes      
c transform t_X, t_Y, and t_Z with r, a, b, and c
      do j = first, last
!if core:
C          if (mod(j-1,n/n_c).eq.0) then
          if ( type(j).ne.0 ) then

             k  = (indexj(j)-1)*t_n/n_c
             m1 = 3*(j-1) + 1
             m2 = m1 + 1
             m3 = m2 + 1

!move its tails
cccccccccccccc
ccc Modification: DURBA and TONI
c     (We clarify the logic of the operation below)
c     The tails will rotate solidary to the core, this means that the
c     dot product of a,b,c and t-r will be invariant under the rotation.
c     This leads to the formula 
c     t_n = r_n + (A_n^-1)A(t-r)
c     where A is a matrix that contains a,b,c on each row, and A_n^-1 is
c     the inverse of the equivalent matrix for a_n,b_n,c_n
c     Since {a,b,c} is orthonormal the matrix A is orthogonal meaning that
c     the inverse and transpose are equivalent, which is very convenient
c     computationally. Of course, this also apply for the new coordinates
c     as far as we ensure the orthonormality of {a_n,b_n,c_n}.
cccccccccccccc
             do i = 1, t_n/n_c
                ki = k + i		
		comp(1) = (t_X(ki)-r(m1))*a(m1) 
     +		   + (t_Y(ki)-r(m2))*a(m2) + (t_Z(ki)-r(m3))*a(m3)		
		comp(2) = (t_X(ki)-r(m1))*b(m1) 
     +		   + (t_Y(ki)-r(m2))*b(m2) + (t_Z(ki)-r(m3))*b(m3)
		comp(3) = (t_X(ki)-r(m1))*c(m1) 
     +		   + (t_Y(ki)-r(m2))*c(m2) + (t_Z(ki)-r(m3))*c(m3)
		   
                t_X_n(ki) = r_n(m1) + a_n(m1)*comp(1)
     +              + b_n(m1)*comp(2) + c_n(m1)*comp(3)
                t_Y_n(ki) = r_n(m2) + a_n(m2)*comp(1)
     +              + b_n(m2)*comp(2) + c_n(m2)*comp(3)
                t_Z_n(ki) = r_n(m3) + a_n(m3)*comp(1)
     +              + b_n(m3)*comp(2) + c_n(m3)*comp(3)
             enddo

ccccccccccc
ccc Modification: TONI (LHref)
ccc Reposition Linker Histones
ccc Since the LH can have flexibility, we need to compute the relative
ccc positions rather than using the static initial positions h_X0
			if (withlink) then
             k  = (indexj(j)-1)*nbLH
             do i = 1,nbLH
                ki = k + i

		comp(1) = (h_X(ki)-r(m1))*a(m1) 
     +		   + (h_Y(ki)-r(m2))*a(m2) + (h_Z(ki)-r(m3))*a(m3)		
		comp(2) = (h_X(ki)-r(m1))*b(m1) 
     +		   + (h_Y(ki)-r(m2))*b(m2) + (h_Z(ki)-r(m3))*b(m3)
		comp(3) = (h_X(ki)-r(m1))*c(m1) 
     +		   + (h_Y(ki)-r(m2))*c(m2) + (h_Z(ki)-r(m3))*c(m3)
		   
                h_X_n(ki) = r_n(m1) + a_n(m1)*comp(1)
     +              + b_n(m1)*comp(2) + c_n(m1)*comp(3)
                h_Y_n(ki) = r_n(m2) + a_n(m2)*comp(1)
     +              + b_n(m2)*comp(2) + c_n(m2)*comp(3)
                h_Z_n(ki) = r_n(m3) + a_n(m3)*comp(1)
     +              + b_n(m3)*comp(2) + c_n(m3)*comp(3)

		
             enddo
			end if



			
          endif
      enddo



      return
      end

