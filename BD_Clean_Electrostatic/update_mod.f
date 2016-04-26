cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Subroutine computes the Euler angles, alpha, beta, gamma, and
c   alpha_p, beta_p, gamma_p, and
c   a_dna, b_dna, c_dna for the chromatin model.
c   length is the array of distances between core/linker bead i and i+1
c
cccccccccccccccccccccccc
ccc Modification: DURBA and TONI
c     In this subroutine a, b, c of linker DNA are changed, but
c     the original formulas have been revised and corrected.
c     In the old code (update subroutine) this was leading to inconsistencies
c     in the MC cycles and errors in the elastic energy.
cccccccccccccccccccccccc

      subroutine update_mod(n_c,nc3,n,n3,type,r,ro,d1,go, a,b,c,
     +                  alpha,beta,gamma, length,a_dna,b_dna,c_dna,
     +                  alpha_p,beta_p,gamma_p,fcount)

C	  use mpi
      implicit NONE

      integer n_c,nc3,n,n3, type(n),myid,ierr,index,j1,j2,j3,fcount
      double precision r(n3),ro,d1,go, a(n3),b(n3),c(n3)
      double precision alpha(n),beta(n),gamma(n), length(n)
      double precision a_dna(nc3),b_dna(nc3),c_dna(nc3)
      double precision alpha_p(n_c),beta_p(n_c),gamma_p(n_c)

      double precision r_forw(3), mi
      double precision da(3), a_old(3)
      double precision a_m(3),b_m(3)
      double precision Ac, apg, f1, f2, ada, bda, si,co
      double precision sa,ca, sb,cb, sg,cg
      double precision R21,R22,R23,R31,R32,R33
      integer i, i1,i2,i3, if1,if2,if3, nm1, ic,ic1,ic2,ic3
      integer nt

      double precision dra,drb,drc
      double precision drperp(3)
      double precision drperp2
      double precision b_old(3),c_old(3)
      double precision a_new(3),b_new(3),c_new(3)
      double precision b_tmp(3)
      double precision prjba

C	  write(*,*) 'a at beginning of update_mod=',myid,a(1)

      character(len=1024) :: filename
      write (filename, "(A20,I3,A5)") "Coordinates data", fcount,'.txt'

       open(unit=1360,name=filename,access='SEQUENTIAL',
     +        status='unknown')
       
C      write (*,*) 'ENTERING UPDATE MOD on PROC' 
C      write(*,*) 'proc,r=',myid,r(1)
C	  write(*,*) 'a=',myid,a(1)
C	  write(*,*) 'b=',myid,b(1)
C	  write(*,*) 'c=',myid,c(1)
      si = dsin(go)
      co = dcos(go)
      nm1 = n-1
      nt = n/n_c

c     This is intended to compute the length (streching) between consecutive beads
c        (i) linker DNA / linker DNA
c        (ii) linker DNA / core
c        (iii) core / linker DNA
C      if (myid.eq.0) then
      do 10 i = 1,nm1
          
        i1 = 3*(i-1) + 1
        i2 = i1 + 1
        i3 = i2 + 1
        
        if1 = i1 + 3
        if2 = if1 + 1
        if3 = if2 + 1
      
        if ( type(i) .eq. 0 ) then

          a_old(1) = a(i1)
          a_old(2) = a(i2)
          a_old(3) = a(i3)
          if ( type(i+1) .eq. 0 ) then
            r_forw(1) = r(if1) - r(i1)
            r_forw(2) = r(if2) - r(i2)
            r_forw(3) = r(if3) - r(i3)
          else
            b_m(1) = -si*a(if1) + co*b(if1)
            b_m(2) = -si*a(if2) + co*b(if2)
            b_m(3) = -si*a(if3) + co*b(if3)
            r_forw(1) = r(if1) - ro*b_m(1) + d1*c(if1) - r(i1)
            r_forw(2) = r(if2) - ro*b_m(2) + d1*c(if2) - r(i2)
            r_forw(3) = r(if3) - ro*b_m(3) + d1*c(if3) - r(i3)
          end if
          length(i) = dsqrt( r_forw(1)**2 +
     +                       r_forw(2)**2 +
     +                       r_forw(3)**2 )
          mi = 1.d0 / length(i)


c     This portion ensures that the a vector of a linker DNA bead points
c     to the position of the next linker DNA bead (when those are contiguous).
c     The b vector changes accordingly based on the previous position
c     (this operation was incorrect in the original code)
c     A Gram-Schmidt process is applied finally to b and c to ensure the orthogonality
c     and normalization of the vectors.

          b_old =(/b(i1),b(i2),b(i3)/)
          c_old =(/c(i1),c(i2),c(i3)/)
          
c     The new a vector (normalized) is related to the distances between consecutives beads
          a_new = mi*r_forw

c     If a vector has changed, then the b vector has to change accordingly
c     In this case we assume that {a,b,c} are orthonormal, and,
c     following this logic, we develop an algorithm that properly reorients b
c     This was based on the algorithm of the old code, but has been corrected and generalized
          if(ALL(a_new.EQ.a_old)) then
             b_new = b_old
             c_new = c_old
          else

c     We compute the projections of r_forw in the old {a,b,c} basis
             dra = dot_product(a_old,r_forw)
             drb = dot_product(b_old,r_forw)
             drc = dot_product(c_old,r_forw)

c     This is the square length of the perpendicular component of r_forw with b
             drperp2 = dra*dra + drc*drc
          
c     If the new a vector is parallel to b
c     (we assume that the two linker DNA beads are not in the same position, i.e., length>0
             if(drperp2.EQ.0) then
                if(drb.GT.0) then
                   b_new = - a_old
                else 
                   b_new = a_old
                end if
             else
c     If the new a vector is oriented in the same hemisphere of the old a vector,
c     then the formula is b = b - (drb/drperp2)*drperp
c     Otherwise is b = - (b - (drb/drperp2)*drperp)
c     This operation in general leads to a non-unitarian b vector.
c     Watch out!!! Although we have treated the cases with drperp2 = 0 separately,
c     the method could have numerical problems if drperp2 is close to zero.
c     However, this is probably unlikely to happen because it would be associated to very 
c     large translations between two beads in a single move.
                drperp = dra*a_old + drc*c_old
                b_tmp = b_old - drb/drperp2*drperp
                if(dra.GE.0) then
                   b_new = b_tmp
                else
                   b_new = -b_tmp
                end if
             end if

c     We ensure orthonormality of b and c
c     First we apply Gram-Shcmidt (necessary or precise enough????)
             prjba = dot_product(b_new,a_new)
             b_new = b_new - prjba*a_new
             b_new = b_new/(sqrt(dot_product(b_new,b_new)))
             
             c_new(1) = a_new(2)*b_new(3) - a_new(3)*b_new(2)
             c_new(2) = a_new(3)*b_new(1) - a_new(1)*b_new(3)
             c_new(3) = a_new(1)*b_new(2) - a_new(2)*b_new(1)
             c_new = c_new/(sqrt(dot_product(c_new,c_new)))
          end if


C	  write(*,*) 'a before assigning new of update_mod=',myid,a(1)


c     We assign the value of the new vectors
          a(i1) = a_new(1)
          a(i2) = a_new(2)
          a(i3) = a_new(3)

          b(i1) = b_new(1)
          b(i2) = b_new(2)
          b(i3) = b_new(3)

          c(i1) = c_new(1)
          c(i2) = c_new(2)
          c(i3) = c_new(3)
          
          


c     Old code
!          a(i1) = mi * r_forw(1)
!          a(i2) = mi * r_forw(2)
!          a(i3) = mi * r_forw(3)

!         da(1) = a(i1) - a_old(1)
!         da(2) = a(i2) - a_old(2)
!         da(3) = a(i3) - a_old(3)

!         bda = b(i1)*da(1) + b(i2)*da(2) + b(i3)*da(3)
!         b(i1) = b(i1) - bda*a_old(1)
!         b(i2) = b(i2) - bda*a_old(2)
!         b(i3) = b(i3) - bda*a_old(3)

!         bda = b(i1)*a(i1) + b(i2)*a(i2) + b(i3)*a(i3)
!         b(i1) = b(i1) - bda*a(i1)
!         b(i2) = b(i2) - bda*a(i2)
!         b(i3) = b(i3) - bda*a(i3)
!         mi = 1.d0/dsqrt( b(i1)**2 +
!    +                     b(i2)**2 +
!    +                     b(i3)**2 )
!         b(i1) = mi * b(i1)
!         b(i2) = mi * b(i2)
!         b(i3) = mi * b(i3)

!         c(i1) = a(i2)*b(i3) - a(i3)*b(i2)
!         c(i2) = a(i3)*b(i1) - a(i1)*b(i3)
!         c(i3) = a(i1)*b(i2) - a(i2)*b(i1)
cccccccccccccccccccccc
        else
          r_forw(1) = r(if1) - ( r(i1)-ro*b(i1)-d1*c(i1) )
          r_forw(2) = r(if2) - ( r(i2)-ro*b(i2)-d1*c(i2) )
          r_forw(3) = r(if3) - ( r(i3)-ro*b(i3)-d1*c(i3) )
          length(i) = dsqrt( r_forw(1)**2 +
     +                       r_forw(2)**2 +
     +                       r_forw(3)**2 )

        end if


   10 continue
   
C	  endif !........................MP correct a,b,c, LOOP END	


C	  call MPI_BCAST(r, n3, MPI_DOUBLE_PRECISION, 
C     +	    0,MPI_COMM_WORLD,ierr)
C      call MPI_BCAST(a, n3, MPI_DOUBLE_PRECISION, 
C     +	    0,MPI_COMM_WORLD,ierr)
C	  call MPI_BCAST(b, n3, MPI_DOUBLE_PRECISION, 
C     +	    0,MPI_COMM_WORLD,ierr)
C	  call MPI_BCAST(c, n3, MPI_DOUBLE_PRECISION, 
C     +	    0,MPI_COMM_WORLD,ierr)


	  ic=0
      do 20 i = 1,nm1
        i1 = 3*(i-1) + 1
        i2 = i1 + 1
        i3 = i2 + 1
        if1 = i1 + 3
        if2 = if1 + 1
        if3 = if2 + 1

        if ( type(i) .eq. 0 ) then
          if ( type(i+1) .eq. 0 ) then
c Regular DNA segment:
            ada = a(i1)*a(if1)+a(i2)*a(if2)+a(i3)*a(if3)
            if ( ada .gt. 1.d0 ) ada = 1.d0
            if ( ada .lt. -1.d0 ) ada = -1.d0
            beta(i) = dacos( ada )
            sb = dsin(beta(i))
            if (beta(i) .ge. 1.0d-10) then
              f1 = (a(if1)*b(i1)+a(if2)*b(i2)+a(if3)*b(i3))/sb
            else
              f1 = (b(if1)*b(i1)+b(if2)*b(i2)+b(if3)*b(i3))
            end if
            if ( f1 .gt. 1.d0 ) f1 = 1.d0
            if ( f1 .lt. -1.d0 ) f1 = -1.d0
            Ac = dacos( f1 )
            f2 = a(if1)*c(i1)+a(if2)*c(i2)+a(if3)*c(i3)
            if ( f2 .ge. 0.0d0 ) then
              alpha(i) = Ac
            else
              alpha(i) = -Ac
            end if
            f1 = ( b(i1)*b(if1)+b(i2)*b(if2)+b(i3)*b(if3) +
     +             c(i1)*c(if1)+c(i2)*c(if2)+c(i3)*c(if3) ) /
     +           ( 1.d0 + ada )
            if ( f1 .gt. 1.d0 ) f1 = 1.d0
            if ( f1 .lt. -1.d0 ) f1 = -1.d0
            apg = dacos( f1 )
            f2 = ( c(i1)*b(if1)+c(i2)*b(if2)+c(i3)*b(if3) -
     +            (b(i1)*c(if1)+b(i2)*c(if2)+b(i3)*c(if3)) ) /
     +           ( 1.0d0 + ada )
            if ( f2 .ge. 0.d0 ) then
              gamma(i) = apg - alpha(i)
            else
              gamma(i) = -apg - alpha(i)
            end if

          else
c DNA segment with bead(i+1) is a core:
            a_m(1) = co*a(if1) + si*b(if1)
            a_m(2) = co*a(if2) + si*b(if2)
            a_m(3) = co*a(if3) + si*b(if3)
            b_m(1) = -si*a(if1) + co*b(if1)
            b_m(2) = -si*a(if2) + co*b(if2)
            b_m(3) = -si*a(if3) + co*b(if3)

            ada = a(i1)*a_m(1)+a(i2)*a_m(2)+a(i3)*a_m(3)
            if ( ada .gt. 1.d0 ) ada = 1.d0
            if ( ada .lt. -1.d0 ) ada = -1.d0
            beta(i) = dacos( ada )
            sb = dsin( beta(i) )
            if (beta(i) .ge. 1.0d-10) then
              f1 = (a_m(1)*b(i1)+a_m(2)*b(i2)+a_m(3)*b(i3))/sb
            else
              f1 = (b_m(1)*b(i1)+b_m(2)*b(i2)+b_m(3)*b(i3))
            end if
            if ( f1 .gt. 1.d0 ) f1 = 1.d0
            if ( f1 .lt. -1.d0 ) f1 = -1.d0
            Ac = dacos( f1 )
            f2 = ( a_m(1)*c(i1)+a_m(2)*c(i2)+a_m(3)*c(i3) )
            if ( f2 .ge. 0.0d0 ) then
              alpha(i) = Ac
            else
              alpha(i) = -Ac
            end if
            f1 = ( b(i1)*b_m(1)+b(i2)*b_m(2)+b(i3)*b_m(3) +
     +             c(i1)*c(if1)+c(i2)*c(if2)+c(i3)*c(if3) ) /
     +           ( 1.d0 + ada )
            if ( f1 .gt. 1.d0 ) f1 = 1.d0
            if ( f1 .lt. -1.d0 ) f1 = -1.d0
            apg = dacos( f1 )
            f2 = ( c(i1)*b_m(1)+c(i2)*b_m(2)+c(i3)*b_m(3) -
     +            (b(i1)*c(if1)+b(i2)*c(if2)+b(i3)*c(if3)) ) /
     +           ( 1.0d0 + ada )
            if ( f2 .ge. 0.d0 ) then
              gamma(i) = apg - alpha(i)
            else
              gamma(i) = -apg - alpha(i)
            end if
c ********  Abotaleb
      write(unit=1360, fmt=9089)'for Linker DNA followed by Core:'
              
          write(unit=1360+fcount, fmt=9809)" am_x(",i,") =",a_m(1),
     +    " am_y(",i,") =",  a_m(2)," am_z(",i,") =", a_m(3)
        write(unit=1360+fcount, fmt=9809)" bm_x(",i,") =",b_m(1),
     +    " bm_y(",i,") =",  b_m(2)," bm_z(",i,") =", b_m(3)
        write(unit=1360+fcount, fmt=9809)" cm_x(",i,") =",c(i1),
     +    " cm_y(",i,") =",  c(i2)," cm_z(",i,") =", c(i3)
       
        
            
            
          end if

        else
c DNA segment with bead(i) is a core:
          ic = ic +1
          ic1 = 3*(ic-1) + 1
          ic2 = ic1+1
          ic3 = ic2+1

          a_dna(ic1) = ( r(if1) - (r(i1)-ro*b(i1)-d1*c(i1)) ) 
     +                  / length(i)
          a_dna(ic2) = ( r(if2) - (r(i2)-ro*b(i2)-d1*c(i2)) )
     +                  / length(i)
          a_dna(ic3) = ( r(if3) - (r(i3)-ro*b(i3)-d1*c(i3)) )
     +                  / length(i)


c %%%%%%%%%%%%%
          cb = a(i1)*a_dna(ic1)+a(i2)*a_dna(ic2)+a(i3)*a_dna(ic3)
          if ( cb .gt. 1.d0 ) cb = 1.d0
          if ( cb .lt. -1.d0 ) cb = -1.d0
          beta_p(ic) = dacos( cb )
          sb = sin( beta_p(ic) )
c
          if ( beta_p(ic) .ge. 1.0d-10 ) then
            b_m(1) = ( a_dna(ic1) - cb*a(i1) ) / sb
            b_m(2) = ( a_dna(ic2) - cb*a(i2) ) / sb
            b_m(3) = ( a_dna(ic3) - cb*a(i3) ) / sb
            ca = b_m(1)*b(i1)+b_m(2)*b(i2)+b_m(3)*b(i3)
            if ( ca .gt. 1.d0 ) ca = 1.d0
            if ( ca .lt. -1.d0 ) ca = -1.d0
            Ac = dacos( ca )
            f1 = a_dna(ic1)*c(i1)+a_dna(ic2)*c(i2)+a_dna(ic3)*c(i3) 
            if (f1 .ge. 0.0d0) then
              alpha_p(ic) = Ac
            else
              alpha_p(ic) = -Ac
            end if
            gamma_p(ic) = -alpha_p(ic)
            sa = dsin(alpha_p(ic))
            sg = dsin(gamma_p(ic))
            cg = dcos(gamma_p(ic))
            R21 = -cg*sb
            R22 = cg*cb*ca-sg*sa
            R23 = cg*cb*sa+sg*ca

            b_dna(ic1) = R21*a(i1) +
     +                   R22*b(i1) +
     +                   R23*c(i1)
            b_dna(ic2) = R21*a(i2) +
     +                   R22*b(i2) +
     +                   R23*c(i2)
            b_dna(ic3) = R21*a(i3) +
     +                   R22*b(i3) +
     +                   R23*c(i3)

            R31 = sg*sb
            R32 = -sg*cb*ca-cg*sa
            R33 = -sg*cb*sa+cg*ca
            c_dna(ic1) = R31*a(i1) +
     +                   R32*b(i1) +
     +                   R33*c(i1)
            c_dna(ic2) = R31*a(i2) +
     +                   R32*b(i2) +
     +                   R33*c(i2)
            c_dna(ic3) = R31*a(i3) +
     +                   R32*b(i3) +
     +                   R33*c(i3)

          else
            b_dna(ic1) = b(i1)
            b_dna(ic2) = b(i2)
            b_dna(ic3) = b(i3)
            c_dna(ic1) = c(i1)
            c_dna(ic2) = c(i2)
            c_dna(ic3) = c(i3)
          end if
c Now get {alpha,beta,gamma} to transform from
c {a,b,c}_dna to {a,b,c}_{i+1}
          ada = a_dna(ic1)*a(if1)+a_dna(ic2)*a(if2)+a_dna(ic3)*a(if3)
          if ( ada .gt. 1.d0 ) ada = 1.d0
          if ( ada .lt. -1.d0 ) ada = -1.d0
          beta(i) = dacos( ada )
          sb = dsin(beta(i))
          if (beta(i) .ge. 1.0d-10) then
            f1 = (a(if1)*b_dna(ic1)+a(if2)*b_dna(ic2)+
     +            a(if3)*b_dna(ic3))/sb
          else
            f1=(b(if1)*b_dna(ic1)+b(if2)*b_dna(ic2)+b(if3)*b_dna(ic3))
          end if
          if ( f1 .gt. 1.d0 ) f1 = 1.d0
          if ( f1 .lt. -1.d0 ) f1 = -1.d0
          Ac = dacos( f1 )
          f2 = a(if1)*c_dna(ic1)+a(if2)*c_dna(ic2)+a(if3)*c_dna(ic3)
          if ( f2 .ge. 0.0d0 ) then
            alpha(i) = Ac
          else
            alpha(i) = -Ac
          end if
          f1=( b_dna(ic1)*b(if1)+b_dna(ic2)*b(if2)+b_dna(ic3)*b(if3) +
     +         c_dna(ic1)*c(if1)+c_dna(ic2)*c(if2)+c_dna(ic3)*c(if3))/
     +       ( 1.d0 + ada )
          if ( f1 .gt. 1.d0 ) f1 = 1.d0
          if ( f1 .lt. -1.d0 ) f1 = -1.d0
          apg = dacos( f1 )
          f2=( c_dna(ic1)*b(if1)+c_dna(ic2)*b(if2)+c_dna(ic3)*b(if3) -
     +        (b_dna(ic1)*c(if1)+b_dna(ic2)*c(if2)+b_dna(ic3)*c(if3)))/
     +       ( 1.0d0 + ada )
          if ( f2 .ge. 0.d0 ) then
            gamma(i) = apg - alpha(i)
          else
            gamma(i) = -apg - alpha(i)
          end if
c ********  Abotaleb
      write(unit=1360+fcount, fmt=9089)'for Core DNA Beads:'
      
          
       write(unit=1360, fmt=9019)" a_dna_x(",i,") =",a_dna(ic1),
     +    " a_dna_y(",i,") =",a_dna(ic2)," a_dna_z(",i,") =",a_dna(ic3)
        write(unit=1360, fmt=9019)" b_dna_x(",i,") =",b_dna(ic1),
     +    " b_dna_y(",i,") =",b_dna(ic2)," b_dna_z(",i,") =",b_dna(ic3)
        write(unit=1360, fmt=9019)" c_dna_x(",i,") =",c_dna(ic1),
     +" c_dna_y(",i,") =",c_dna(ic2)," c_dna_z(",i,") =",c_dna(ic3)
        write(unit=1360, fmt=9099)" alpha_p(",i,") =",
     +                            alpha_p(ic)
        write(unit=1360, fmt=9099)" beta_p(",i,") =",
     +                              beta_p(ic)
        write(unit=1360, fmt=9099)" gamma_p(",i,") =",
     +                               gamma_p(ic)
       
          
          
          
        end if

   20 continue
   
C      write (*,*) 'AFTER UPDATE MOD'
C      write(*,*) 'proc,r=',myid,r(1)
C	  write(*,*) 'a=',myid,a(1)
C	  write(*,*) 'b=',myid,b(1)
C	  write(*,*) 'c=',myid,c(1)
C	  write(*,*) 'a at end of update_mod=',myid,a(1)

c ****************************************************
c *****************************************************
c Abotaleb 26 October 2015
c *****************************************************
c *****************************************************
c *****************************************************

       index=0
       write(unit=1360+fcount, fmt=9089)'DNA Beads Data:'
      do 2230 i = 1, n
         j1 = 3*index + 1
         j2 = 3*index + 2
         j3 = 3*index + 3
        write(unit=1360, fmt=9009)" rx(",index+1,") =",r(j1),
     +    " ry(",index+1,") =",  r(j2)," rz(",index+1,") =", r(j3)
        write(unit=1360, fmt=9009)" ax(",index+1,") =",a(j1),
     +    " ay(",index+1,") =",  a(j2)," az(",index+1,") =", a(j3)
        write(unit=1360, fmt=9009)" bx(",index+1,") =",b(j1),
     +    " by(",index+1,") =",  b(j2)," bz(",index+1,") =", b(j3)
        write(unit=1360, fmt=9009)" cx(",index+1,") =",c(j1),
     +    " cy(",index+1,") =",  c(j2)," cz(",index+1,") =", c(j3)
        write(unit=1360,fmt=9099)" alpha(",index+1,") =",
     +    alpha(index+1)
        write(unit=1360, fmt=9099)" beta(",index+1,") =",
     +    beta(index+1)
        write(unit=1360, fmt=9099)" gamma(",index+1,") =",
     + gamma(index+1)
         index=index+1
 2230 continue

      
        close(unit=1360+fcount)
      return
 9809 format(1x, A7, I2, 2x, A3,2X,F10.6,
     +       1x, A7, I2, 2x, A3,2x,F10.6,
     +       1x, A7, I2, 2x, A3,2x,F14.7)
 9009 format(1x, A4, I2, 2x, A3,2X,F10.6,
     +       1x, A4, I2, 2x, A3,2x,F10.6,
     +       1x, A4, I2, 2x, A3,2x,F14.7)
9099  format(1x,A12,I3,2x,A4,2x,F10.6)
 9019 format(1x, A10, I2, 2x, A3,2X,F10.6,
     +       1x, A10, I2, 2x, A3,2x,F10.6,
     +       1x, A10, I2, 2x, A3,2x,F14.7)
     
9089  format(1x,A45)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

