Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
      subroutine geteuler(n_c,nc3,n,n3,type,r,ro,d1,go, a,b,c,
     +                  alpha,beta,gamma, length,a_dna,b_dna,c_dna,
     +                  alpha_p,beta_p,gamma_p)

      implicit none  !CSUN 6/2/04
      
      integer n_c,nc3,n,n3, type(n)
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
      integer nbead

      si = dsin(go)
      co = dcos(go)
      nm1 = n-1
      nt = 6
      nbead=n/n_c

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

        else
          r_forw(1) = r(if1) - ( r(i1)-ro*b(i1)-d1*c(i1) )
          r_forw(2) = r(if2) - ( r(i2)-ro*b(i2)-d1*c(i2) )
          r_forw(3) = r(if3) - ( r(i3)-ro*b(i3)-d1*c(i3) )
          length(i) = dsqrt( r_forw(1)**2 +
     +                       r_forw(2)**2 +
     +                       r_forw(3)**2 )

        end if
   10 continue

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
            if ( f2 .ge. 0. ) then
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
            if ( f2 .ge. 0. ) then
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
          end if

        else
c DNA segment with bead(i) is a core:
COLD 5/14/04          ic = (i-1)/nt + 1
          ic = (i-1)/nbead + 1
          ic1 = 3*(ic-1) + 1
          ic2 = ic1+1
          ic3 = ic2+1

          a_dna(ic1) = ( r(if1) - (r(i1)-ro*b(i1)-d1*c(i1)) ) 
     +                  / length(i)
          a_dna(ic2) = ( r(if2) - (r(i2)-ro*b(i2)-d1*c(i2)) )
     +                  / length(i)
          a_dna(ic3) = ( r(if3) - (r(i3)-ro*b(i3)-d1*c(i3)) )
     +                  / length(i)

          cb = a(i1)*a_dna(ic1)+a(i2)*a_dna(ic2)+a(i3)*a_dna(ic3)
          if ( cb .gt. 1.d0 ) cb = 1.d0
          if ( cb .lt. -1.d0 ) cb = -1.d0
          beta_p(ic) = dacos( cb )
COLD          sb = sin( beta_p(ic) )
          sb = dsin( beta_p(ic) )   ! CSUN 12/8/03
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
            if (f1 .ge. 0.) then
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
          if ( f2 .ge. 0. ) then
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

        end if

   20 continue

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

