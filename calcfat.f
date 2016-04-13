coccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine calcfat(n_c,nc3,n,n3,type,r,a,b,c,
     +                   alpha,beta,gamma,length,
     +                   a_dna,b_dna,c
     +                   alpha_p,beta_p,gamma_p,
     +                   lo,ro,d1,go,h,g,s,k_e,k_ex,
     +                    Nq,Nq3,core_pos,core_q,
     +                    debye,q_l,phi_o,force,torque)

 
      integer n_c,nc3,n,n3, type(n)
      double precision r(n3), a(n3),b(n3),c(n3)
      double precision alpha(n),beta(n),gamma(n), length(n)
      double precision a_dna(nc3),b_dna(nc3),c_dna(nc3)
      double precision alpha_p(n_c),beta_p(n_c),gamma_p(n_c)
      double precision lo,ro,d1,go,h,g,s,k_e,k_ex
      integer Nq, Nq3
      double precision core_pos(Nq3), core_q(Nq)
      double precision debye, q_l, phi_o
      double precision force(n3), torque(n3)

      integer i,j,k,l, nt
      integer i1,i2,i3, ib1,ib2,ib3, if1,if2,if3
      integer nm1,nm2, im1,im2,im3, ic, ic1,ic2,ic3
      integer j1,j2,j3, k1,k2,k3, l1,l2,l3
      double precision c1,c2, s1,s2, g1,g2, dist
      double precision a_m(3), si,co
      double precision df(3)
      double precision Stri(3), Strim1(3)
      double precision Ai(3), Aim1(3), Bi(3), Bim1(3)
      double precision Chi(3),Chim1(3), Zhi(3),Zhim1(3)
      double precision t_temp(3), mag, fa,fb,fc
      double precision mi, z(3), ql_ql
      double precision ada,adb,adc, cda,bda,cdb,bdb,cdc,bdc
      double precision  Rcut

      nt = 6
      nm1 = n-1
      nm2 = nm1-1
      Rcut = 22.0d0

      si = dsin( go )
      co = dcos( go )     
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Bending on 1st partical [core]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       4 Oct/2015  [Old]


      if ( beta(1) .ge. 1.0d-10 ) then
        g1 = beta(1) / ( dsin( beta(1) )*length(1) )
      else
        g1 = 1.0d0 / length(1)
      end if
      
      c1 = dcos( beta(1) )
      Ai(1) = g1*( a(4)-c1*a_dna(1) )
      Ai(2) = g1*( a(5)-c1*a_dna(2) )
      Ai(3) = g1*( a(6)-c1*a_dna(3) )
      
      if ( beta_p(1) .ge. 1.0d-10 ) then
        g2 = beta_p(1) / ( dsin( beta_p(1) )*length(1) )
      else 
        g2 = 1.0d0 / length(1)
      end if
      
      c2 = dcos( beta_p(1) )
      Bi(1) = g2*( a(1)-c2*a_dna(1) )
      Bi(2) = g2*( a(2)-c2*a_dna(2) )
      Bi(3) = g2*( a(3)-c2*a_dna(3) )

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Bending on 1st partical [core]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       4 Oct/2015 [Old]

      g1 = (alpha(1)+gamma(1)-phi_o)*dtan(0.5d0*beta(1)) /
     +      length(1)
      c1 = dcos( alpha(1) )
      s1 = dsin( alpha(1) )
      Chi(1) = g1*( c1*c_dna(1)-s1*b_dna(1) )
      Chi(2) = g1*( c1*c_dna(2)-s1*b_dna(2) )
      Chi(3) = g1*( c1*c_dna(3)-s1*b_dna(3) )

      g2 = (alpha(1)+gamma(1)-phi_o)*dtan(0.5d0*beta_p(1)) / 
     +      length(1)
      c2 = dcos( gamma_p(1) )
      s2 = dsin( gamma_p(1) )
      Zhi(1) = g2*( c2*c_dna(1)+s2*b_dna(1) )
      Zhi(2) = g2*( c2*c_dna(2)+s2*b_dna(2) )
      Zhi(3) = g2*( c2*c_dna(3)+s2*b_dna(3) )

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Stretching on 1st partical [core]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       4 Oct/2015 [Old]

      Stri(1) = (length(1)-lo)*a_dna(1)
      Stri(2) = (length(1)-lo)*a_dna(2)
      Stri(3) = (length(1)-lo)*a_dna(3)

      Strim1(1) = 0.d0
      Strim1(2) = 0.d0
      Strim1(3) = 0.d0


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
CCCCCCCCCCCCC Total Mechanical force on 1st bead [Core]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       4 Oct/2015 [Old]
      
      force(1) = h*( Stri(1) )
     +         - g*( Ai(1) + Bi(1) )
     +         + s*( Chi(1) + Zhi(1) )
      force(2) = h*( Stri(2) )
     +         - g*( Ai(2) + Bi(2) )
     +         + s*( Chi(2) + Zhi(2) )
      force(3) = h*( Stri(3) )
     +         - g*( Ai(3) + Bi(3) )
     +         + s*( Chi(3) + Zhi(3) )
      
      
      
      
      end