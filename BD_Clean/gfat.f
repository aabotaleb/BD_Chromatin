coccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gfat(n_c,nc3,n,n3,type,r,a,b,c,E,
     +                alpha,beta,gamma,length,
     +                a_dna,b_dna,c_dna,
     +                alpha_p,beta_p,gamma_p,
     +                lo,ro,d1,go,h,g,s,k_e,k_ex,hd2,gd2,sd2,
     +                Nq,Nq3,core_pos,core_q,
     +                debye,q_l,phi_o,force,torque,t_n,t_n3, ! t_n,t_n3 number of tails and coors
     +                force_dnaStr_X,force_dnaStr_Y,force_dnaStr_Z,
     +                force_dnaBend_X,force_dnaBend_Y,force_dnaBend_Z,
     +             force_dnaTwist_X,force_dnaTwist_Y,force_dnaTwist_Z,
     +            h_X,h_Y,h_Z,h_chg,h_n,
     +                t_grp,t_fix,t_chg,t_rad,t_mass,        ! tail business
     +                t_X,t_Y,t_Z, t_X0,t_Y0,t_Z0,           ! >>
     +                t_bond,t_bond_v,t_bond_c,              ! >>
     +                t_angle, t_angle_v, t_angle_c,force_t, ! >>
     +                force_tStrInt_X,force_tStrInt_Y,force_tStrInt_Z,
     +                force_tBend_X,force_tBend_Y,force_tBend_Z,    ! isolating internal stretching and bending forces
     +                t_Eb_num,t_Ea_num,Eb_tc_num,           ! Numerical values of stretching,bending energies in tail beads to double check
     +                 np,myid,ierr)             

c Adding all tail parameters 
      
      integer n_c,nc3,n,n3, type(n)
      integer t_n,t_n3    ! Abotaleb
      double precision r(n3), a(n3),b(n3),c(n3)
      double precision alpha(n),beta(n),gamma(n), length(n)
      double precision a_dna(nc3),b_dna(nc3),c_dna(nc3)
      double precision alpha_p(n_c),beta_p(n_c),gamma_p(n_c)
      double precision lo,ro,d1,go,h,g,s,k_e,k_ex
      integer Nq, Nq3
      double precision core_pos(Nq3), core_q(Nq)
      double precision debye, q_l, phi_o(n_c)
      double precision force(n3), torque(n3)
      ! Abotaleb  added to account for the forces due to tail beads
      double precision force_t(t_n3),force_tStrInt_X(t_n)
      double precision force_tStrInt_Y(t_n),force_tStrInt_Z(t_n)
      double precision force_tBend_X(t_n),force_tBend_Y(t_n)
      double precision force_tBend_Z(t_n)
      double precision t_Eb_num,t_Ea_num,Eb_tc_num
      
      integer i,j,k,l, nt,coreInd
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
ccccccccc tail params    - Abotaleb
      double precision t_chg(t_n), t_rad(t_n), t_mass(t_n) 
      double precision t_X(t_n), t_Y(t_n), t_Z(t_n)
      double precision t_X0(t_n), t_Y0(t_n), t_Z0(t_n)
      double precision t_bond(t_n), t_bond_v(t_n), t_bond_c(t_n)
      double precision t_angle(t_n), t_angle_v(t_n), t_angle_c(t_n)
      integer t_grp(t_n), t_fix(t_n)
      double precision dx, dy, dz, dr2, dr, dt2
      double precision x12, y12, z12, r12, r12s
      double precision x10, y10, z10, r10, r10s
      double precision phi
      double precision temp
      
      !Abotaleb added to account for isolated force X,Y,Z and for each force type 
      double precision force_dnaStr_X(n),force_dnaStr_Y(n)
      double precision force_dnaStr_Z(n)
      double precision force_dnaBend_X(n),force_dnaBend_Y(n)
      double precision force_dnaBend_Z(n)
      double precision force_dnaTwist_X(n),force_dnaTwist_Y(n)
      double precision force_dnaTwist_Z(n)
c     local variables needed for testgh part
      double precision fnc, yhy, xc(n3),gc(n3),y(n3),vec(n3)
      integer isAnalytic      
      double precision E(6)
      double precision aold(n3),bold(n3),cold(n3)
      double precision hd2,gd2,sd2
cccccccc paralleization variables - Abotaleb [tail params only added]
        integer np,myid,ierr,t_interval
        integer t_par_start,t_par_stop
c        Variables needed for the potential calculation routine
      
       
      double precision Ev,Ec , ssle,ssbe,sstw
      double precision   Rcut2
      integer h_n
      double precision   h_chg(h_n)
       

      double precision h_X(h_n), h_Y(h_n), h_Z(h_n)
     
       
      double precision t_exv_e,t_exv_d,evd_hcG,evd_hcC,evd_hlG,evd_hlC
      double precision t_Eb, t_Ea, t_Ec, t_Ev
      
      double precision vdw_cut, vdw_cut2
      double precision evd_tc, evd_tl, evd_cc, evd_cl, evd_ll, evd_hh
      integer m1, m2, m3

      double precision c0 
      double precision test1, test2, test3
      double precision SMALL, SMALL2
      parameter(SMALL = 1.0d-4, SMALL2 = 1.0d-8)
c ... The 6 components of electrostatics, L: linker; C: core; T: tail
      double precision Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL, h_tc
c ... The 2 subcomponents of Ec_TT: Ec_TT1: tails in the same core; Ec_TT2: tails in different cores
c     The 2 subcomponents of Ec_CT: Ec_CT1: tails with its own core; Ec_CT2: tails with other cores 
      double precision Ec_TT1, Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH
      integer ntpc
      double precision Eb_tc,x0,y0,z0,evd_link,Ec_HL1,Ec_HL2,Ec_HC
      logical	withlink
      double precision  spring1_E, spring2_E
      double precision  fs_x, fs_y, fs_z, ss_x, ss_y, ss_z
      double precision  dist_s1, dist_s2
      double precision  first_s, second_s


      integer fnonparLH
      double precision debyell
      double precision Ev_LL, Ev_CC, Ev_TT, Ev_LC, Ev_CT, Ev_TL
      double precision Ev_TT1, Ev_TT2, Ev_CT1, Ev_CT2
      double precision Ev_HH, Ev_HT, Ev_HL1, Ev_HL2, Ev_HC
      double precision e_tmp,e_tmp1,e_tmp2
      common /coeff/ debyell,fnonparLH !!! Common block for coefficients
      common /Ev/ Ev_LL,Ev_CC,Ev_TT,Ev_LC,Ev_CT,Ev_TL,Ev_TT1,Ev_TT2,
     +     Ev_CT1,Ev_CT2,Ev_HH,Ev_HT,Ev_HL1,Ev_HL2,Ev_HC

      integer jlow,jhigh,grpi,grpj
      double precision rehh,rehc,rehl,revd
      double precision kstr,kben,ktor
      double precision kstrd2,kbend2,ktord2,ktmp
      double precision ELHel
      !double precision ELHstr,ELHben,ELHtor (global)
      integer k0
      integer nconn,nbend
      
      integer indexj(n)
          open(unit=199,name='Analytic Forces[MATLAB Version].txt',
     +access='SEQUENTIAL',        status='unknown')
        
      
      nt = 6
      nm1 = n-1
      nm2 = nm1-1
      Rcut = 22.0d0

      si = dsin( go )
      co = dcos( go )

      
cccccccc       1st     Part                      cccccccc
      
c First, is calculated the mechanical forces and torques:
c  Stretching, Bending, Twisting forces and the associated
c  torques. After this, the bead-bead (N^2) forces and
c  torques are calculated.
c
c Bending on first paricle (a core):
c
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
c
c Twisting force on first particle (a core)
c
      g1 = (alpha(1)+gamma(1)-phi_o(1))*dtan(0.5d0*beta(1)) /
     +      length(1)
      c1 = dcos( alpha(1) )
      s1 = dsin( alpha(1) )
      Chi(1) = g1*( c1*c_dna(1)-s1*b_dna(1) )
      Chi(2) = g1*( c1*c_dna(2)-s1*b_dna(2) )
      Chi(3) = g1*( c1*c_dna(3)-s1*b_dna(3) )

      g2 = (alpha(1)+gamma(1)-phi_o(1))*dtan(0.5d0*beta_p(1)) / 
     +      length(1)
      c2 = dcos( gamma_p(1) )
      s2 = dsin( gamma_p(1) )
      Zhi(1) = g2*( c2*c_dna(1)+s2*b_dna(1) )
      Zhi(2) = g2*( c2*c_dna(2)+s2*b_dna(2) )
      Zhi(3) = g2*( c2*c_dna(3)+s2*b_dna(3) )
c
c Stretching force on first particle (a core)
c
      Stri(1) = (length(1)-lo)*a_dna(1)
      Stri(2) = (length(1)-lo)*a_dna(2)
      Stri(3) = (length(1)-lo)*a_dna(3)

      Strim1(1) = 0.d0
      Strim1(2) = 0.d0
      Strim1(3) = 0.d0
c
      force(1) = h*( Stri(1) )
     +         - g*( Ai(1) + Bi(1) )
     +         + s*( Chi(1) + Zhi(1) )
      force(2) = h*( Stri(2) )
     +         - g*( Ai(2) + Bi(2) )
     +         + s*( Chi(2) + Zhi(2) )
      force(3) = h*( Stri(3) )
     +         - g*( Ai(3) + Bi(3) )
     +         + s*( Chi(3) + Zhi(3) )
c  torques due to forces:
      fa = a(1)*force(1) + a(2)*force(2) + a(3)*force(3)
      fb = b(1)*force(1) + b(2)*force(2) + b(3)*force(3)
      fc = c(1)*force(1) + c(2)*force(2) + c(3)*force(3)
      torque(1) = -ro*fc + d1*fb
      torque(2) = -d1*fa
      torque(3) =  ro*fa
cccccccc isolated force for the first core [first bead is always a core]:
                      force_dnaStr_X(1)=h*( Stri(1) )
                      force_dnaStr_Y(1)=h*( Stri(2) )
                      force_dnaStr_Z(1)=h*( Stri(3) )
                      force_dnaBend_X(1)=- g*( Ai(1) + Bi(1) )
                      force_dnaBend_Y(1)=- g*( Ai(2) + Bi(2) )
                      force_dnaBend_Z(1)=- g*( Ai(3) + Bi(3) )
                      force_dnaTwist_X(1)= s*( Chi(1) + Zhi(1) )
                      force_dnaTwist_Y(1)= s*( Chi(2) + Zhi(2) )
                      force_dnaTwist_Z(1)= s*( Chi(3) + Zhi(3) )

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C             Iterating over the rest beads [Mechanical forces]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      
      coreInd=1;
      do 100 i = 2,nm1
        im1 = i-1
        i1 = 3*im1 + 1
        i2 = i1 + 1
        i3 = i2 + 1
        ib1 = i1 - 3
        ib2 = ib1 + 1
        ib3 = ib2 + 1
        if1 = i1 + 3
        if2 = if1 + 1
        if3 = if2 + 1
c
c (1.) Bending business:
c
        Aim1(1) = Ai(1)
        Aim1(2) = Ai(2)
        Aim1(3) = Ai(3)
        Bim1(1) = Bi(1)
        Bim1(2) = Bi(2)
        Bim1(3) = Bi(3)

c
        if ( type(i) .eq. 0 ) then

c  i-bead is a Linker Bead:
          if ( type(i+1) .eq. 0 ) then
C             and the following is also linker bead
              if ( beta(i) .ge. 1.0d-10 ) then
                 g1 = beta(i) / ( dsin( beta(i) )*length(i) )
              else 
                 g1 = 1.0d0 / length(i)
              end if
            
              c1 = dcos( beta(i) )
              Ai(1) = g1*( a(if1)-c1*a(i1) )
              Ai(2) = g1*( a(if2)-c1*a(i2) )
              Ai(3) = g1*( a(if3)-c1*a(i3) )
              
          else
           ! linker followed by a core
            a_m(1) =  co*a(if1) + si*b(if1)
            a_m(2) =  co*a(if2) + si*b(if2)
            a_m(3) =  co*a(if3) + si*b(if3)
            if ( beta(i) .ge. 1.0d-10 ) then
              g1 = beta(i) / ( dsin( beta(i) )*length(i) )
            else
              g1 = 1.0d0 / length(i)
            end if
            c1 = dcos( beta(i) )
            Ai(1) = g1*( a_m(1)-c1*a(i1) )
            Ai(2) = g1*( a_m(2)-c1*a(i2) )
            Ai(3) = g1*( a_m(3)-c1*a(i3) )
          end if

          if ( type(i-1) .eq. 0 ) then
              ! previous is linker
            if ( beta(im1) .ge. 1.0d-10 ) then
              g2 = beta(im1) / ( dsin( beta(im1) )*length(i) )
            else
              g2 = 1.0d0 / length(i)
            end if
            
            c2 = dcos( beta(im1) )
            Bi(1) = g2*( a(ib1)-c2*a(i1) )
            Bi(2) = g2*( a(ib2)-c2*a(i2) )
            Bi(3) = g2*( a(ib3)-c2*a(i3) )
          else
          ! previous is core        
            ic = (i-2)/nt + 1
            ic1 = 3*(ic-1)+1
            ic2 = ic1+1
            ic3 = ic2+1
            if ( beta(im1) .ge. 1.0d-10 ) then
              g2 = beta(im1) / ( dsin( beta(im1) )*length(i) )
            else
              g2 = 1.0d0 / length(i)
            end if
            c2 = dcos( beta(im1) )
            Bi(1) = g2*( a_dna(ic1)-c2*a(i1) )
            Bi(2) = g2*( a_dna(ic2)-c2*a(i2) )
            Bi(3) = g2*( a_dna(ic3)-c2*a(i3) )
          end if
        else
c
c  i-bead is a Core Bead:
c
          ic = (i-1)/nt + 1
          ic1 = 3*(ic-1)+1
          ic2 = ic1+1
          ic3 = ic2+1

          if ( beta(i) .ge. 1.0d-10 ) then
            g1 = beta(i) / ( dsin( beta(i) )*length(i) )
          else
            g1 = 1.0d0 / length(i)
          end if
          c1 = dcos( beta(i) )
          Ai(1) = g1*( a(if1)-c1*a_dna(ic1) )
          Ai(2) = g1*( a(if2)-c1*a_dna(ic2) )
          Ai(3) = g1*( a(if3)-c1*a_dna(ic3) )
          if ( beta_p(ic) .ge. 1.0d-10 ) then
            g2 = beta_p(ic) / ( dsin( beta_p(ic) )*length(i) )
          else
            g2 = 1.0d0 / length(i)
          end if
          c2 = dcos( beta_p(ic) )
          Bi(1) = g2*( a(i1)-c2*a_dna(ic1) )
          Bi(2) = g2*( a(i2)-c2*a_dna(ic2) )
          Bi(3) = g2*( a(i3)-c2*a_dna(ic3) )

       end if
c  Now Ai(1),Ai(2),Ai(3)  and Bi(1),Bi(2)m,i(3) are calculated for all cases

c
c (2.) Twisting business:
c
        Chim1(1) = Chi(1)
        Chim1(2) = Chi(2)
        Chim1(3) = Chi(3)
        Zhim1(1) = Zhi(1)
        Zhim1(2) = Zhi(2)
        Zhim1(3) = Zhi(3)

        g1 = (alpha(i)+gamma(i)-phi_o(coreInd))*dtan(0.5d0*beta(i)) /
     +        length(i)
        c1 = dcos( alpha(i) )
        s1 = dsin( alpha(i) )
        if ( type(i) .eq. 0 ) then
          g2 = (alpha(im1)+gamma(im1)-phi_o(coreInd))*
     +          dtan(0.5d0*beta(im1)) /length(i)
          
          c2 = dcos( gamma(im1) )
          s2 = dsin( gamma(im1) )
          Chi(1) = g1*( c1*c(i1)-s1*b(i1) )
          Chi(2) = g1*( c1*c(i2)-s1*b(i2) )
          Chi(3) = g1*( c1*c(i3)-s1*b(i3) )
          Zhi(1) = g2*( c2*c(i1)+s2*b(i1) )
          Zhi(2) = g2*( c2*c(i2)+s2*b(i2) )
          Zhi(3) = g2*( c2*c(i3)+s2*b(i3) )
        else
          ic = (i-1)/nt + 1
          ic1 = 3*(ic-1)+1
          ic2 = ic1+1
          ic3 = ic2+1
          g2 = (alpha(i)+gamma(i)-phi_o(coreInd))
     +         *dtan(0.5d0*beta_p(ic)) /length(i)
          c2 = dcos( gamma_p(ic) )
          s2 = dsin( gamma_p(ic) )
          Chi(1) = g1*( c1*c_dna(ic1)-s1*b_dna(ic1) )
          Chi(2) = g1*( c1*c_dna(ic2)-s1*b_dna(ic2) )
          Chi(3) = g1*( c1*c_dna(ic3)-s1*b_dna(ic3) )
          Zhi(1) = g2*( c2*c_dna(ic1)+s2*b_dna(ic1) )
          Zhi(2) = g2*( c2*c_dna(ic2)+s2*b_dna(ic2) )
          Zhi(3) = g2*( c2*c_dna(ic3)+s2*b_dna(ic3) )
        end if
c
c (3.) Stretching business
c
        Strim1(1) = Stri(1)
        Strim1(2) = Stri(2)
        Strim1(3) = Stri(3)

        if ( type(i) .eq. 0 ) then
          Stri(1) = (length(i)-lo)*a(i1)
          Stri(2) = (length(i)-lo)*a(i2)
          Stri(3) = (length(i)-lo)*a(i3)
        else
          ic = (i-1)/nt + 1
          ic1 = 3*(ic-1)+1
          ic2 = ic1+1
          ic3 = ic2+1
          Stri(1) = (length(i)-lo)*a_dna(ic1)
          Stri(2) = (length(i)-lo)*a_dna(ic2)
          Stri(3) = (length(i)-lo)*a_dna(ic3)
        end if
c
c Twisting, bending, stretching forces:
c
c forward segment...
        df(1) = h*( Stri(1) )
     +            - g*( Ai(1)+Bi(1) )
     +            + s*( Chi(1)+Zhi(1) )
        df(2) = h*( Stri(2) )
     +            - g*( Ai(2)+Bi(2) )
     +            + s*( Chi(2)+Zhi(2) )
        df(3) = h*( Stri(3) )
     +            - g*( Ai(3)+Bi(3) )
     +            + s*( Chi(3)+Zhi(3) )
        force(i1) = df(1)
        force(i2) = df(2)
        force(i3) = df(3)
        if ( type(i) .ne. 0 ) then
c torques due to forces
          fa = a(i1)*df(1) + a(i2)*df(2) + a(i3)*df(3)
          fb = b(i1)*df(1) + b(i2)*df(2) + b(i3)*df(3)
          fc = c(i1)*df(1) + c(i2)*df(2) + c(i3)*df(3)
          torque(i1) = -ro*fc + d1*fb
          torque(i2) = -d1*fa
          torque(i3) =  ro*fa
        end if
c back segment...
        df(1) = h*(  -Strim1(1) )
     +            - g*( -Aim1(1)-Bim1(1) )
     +            + s*( -Chim1(1)-Zhim1(1) )
        df(2) = h*( -Strim1(2) )
     +            - g*( -Aim1(2)-Bim1(2) )
     +            + s*( -Chim1(2)-Zhim1(2) )
        df(3) = h*( -Strim1(3) )
     +            - g*( -Aim1(3)-Bim1(3) )
     +            + s*( -Chim1(3)-Zhim1(3) )
        force(i1) = force(i1) + df(1)
        force(i2) = force(i2) + df(2)
        force(i3) = force(i3) + df(3)
        if ( type(i) .ne. 0 ) then
c torques due to forces
          fa = a(i1)*df(1) + a(i2)*df(2) + a(i3)*df(3)
          fb = b(i1)*df(1) + b(i2)*df(2) + b(i3)*df(3)
          fc = c(i1)*df(1) + c(i2)*df(2) + c(i3)*df(3)
          torque(i1) = torque(i1) - co*ro*fc - d1*fb
          torque(i2) = torque(i2) - si*ro*fc + d1*fa
          torque(i3) = torque(i3) + si*ro*fb + co*ro*fa
       end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  

                      force_dnaStr_X(i)=h*( Stri(1) )+h*(  -Strim1(1) )
                      force_dnaStr_Y(i)=h*( Stri(2) )+h*(  -Strim1(2) )
                      force_dnaStr_Z(i)=h*( Stri(3) )+h*(  -Strim1(3) )
      force_dnaBend_X(i)=- g*( Ai(1) + Bi(1) )- g*( -Aim1(1)-Bim1(1) )
      force_dnaBend_Y(i)=- g*( Ai(2) + Bi(2) )- g*( -Aim1(2)-Bim1(2) )
      force_dnaBend_Z(i)=- g*( Ai(3) + Bi(3) )- g*( -Aim1(3)-Bim1(3) )
      force_dnaTwist_X(i)=s*(Chi(1)+Zhi(1)) + s*( -Chim1(1)-Zhim1(1) )
      force_dnaTwist_Y(i)=s*(Chi(2)+Zhi(2)) + s*( -Chim1(2)-Zhim1(2) )
      force_dnaTwist_Z(i)=s*(Chi(3)+Zhi(3)) + s*( -Chim1(3)-Zhim1(3) )

  100 continue
      i1 = n3-2
      i2 = i1+1
      i3 = i2+1
c
c Bending, twisting, stretching forces on bead-n (a dna bead):
c
      force(i1) = -h*( Stri(1) )
     +          + g*( Ai(1) + Bi(1) )
     +          - s*( Chi(1) + Zhi(1) )
      force(i2) = -h*( Stri(2) )
     +          + g*( Ai(2) + Bi(2) )
     +          - s*( Chi(2) + Zhi(2) )
      force(i3) = -h*( Stri(3) )
     +          + g*( Ai(3) + Bi(3) )
     +          - s*( Chi(3) + Zhi(3) )

                      force_dnaStr_X(n)=h*( Stri(1) ) 
                      force_dnaStr_Y(n)=h*( Stri(2) )
                      force_dnaStr_Z(n)=h*( Stri(3) ) 
      force_dnaBend_X(n)=- g*( Ai(1) + Bi(1) )
      force_dnaBend_Y(n)=- g*( Ai(2) + Bi(2) )
      force_dnaBend_Z(n)=- g*( Ai(3) + Bi(3) )
      force_dnaTwist_X(n)=s*(Chi(1)+Zhi(1))
      force_dnaTwist_Y(n)=s*(Chi(2)+Zhi(2))
      force_dnaTwist_Z(n)=s*(Chi(3)+Zhi(3)) 
      
c
c Now, calculate the mechanical Torques:
c
      ada = a_dna(1)*a(1)+a_dna(2)*a(2)+a_dna(3)*a(3)
      adb = a_dna(1)*b(1)+a_dna(2)*b(2)+a_dna(3)*b(3)
      adc = a_dna(1)*c(1)+a_dna(2)*c(2)+a_dna(3)*c(3)

      mag = s*( alpha(1)+gamma(1)-phi_o(coreInd) )
      torque(1) = torque(1) + mag*ada
      torque(2) = torque(2) + mag*adb
      torque(3) = torque(3) + mag*adc
c  extra Bending torques:
      torque(2) = torque(2) - g*beta_p(1)*adc / dsin(beta_p(1))
      torque(3) = torque(3) + g*beta_p(1)*adb / dsin(beta_p(1))
c  extra Twisting torques:
      s1 = dsin( alpha_p(1) )
      c1 = dcos( alpha_p(1) )

      mag = s*( alpha(1)+gamma(1)-phi_o (coreInd))
     +      *dtan( 0.5d0*beta_p(1) )
      cda = c_dna(1)*a(1)+c_dna(2)*a(2)+c_dna(3)*a(3)
      bda = b_dna(1)*a(1)+b_dna(2)*a(2)+b_dna(3)*a(3)
      cdb = c_dna(1)*b(1)+c_dna(2)*b(2)+c_dna(3)*b(3)
      bdb = b_dna(1)*b(1)+b_dna(2)*b(2)+b_dna(3)*b(3)
      cdc = c_dna(1)*c(1)+c_dna(2)*c(2)+c_dna(3)*c(3)
      bdc = b_dna(1)*c(1)+b_dna(2)*c(2)+b_dna(3)*c(3)
      torque(1) = torque(1) - mag*(s1*cda + c1*bda)
      torque(2) = torque(2) - mag*(s1*cdb + c1*bdb)
      torque(3) = torque(3) - mag*(s1*cdc + c1*bdc)

      do 200 i = 2,nm1
        im1 = i-1
        i1 = 3*im1 + 1
        i2 = i1 + 1
        i3 = i2 + 1
        ib1 = i1-3
        ib2 = ib1+1
        ib3 = ib2+1
        if ( type(i) .eq. 0 ) then
          torque(i1) = s*( alpha(i)+gamma(i) -
     +                     alpha(im1)-gamma(im1) )
          torque(i2) = 0.d0
          torque(i3) = 0.d0
        else
          ic = (i-1)/nt + 1
          ic1 = 3*(ic-1)+1
          ic2 = ic1+1
          ic3 = ic2+1

          ada = a_dna(ic1)*a(i1)+a_dna(ic2)*a(i2)+a_dna(ic3)*a(i3)
          adb = a_dna(ic1)*b(i1)+a_dna(ic2)*b(i2)+a_dna(ic3)*b(i3)
          adc = a_dna(ic1)*c(i1)+a_dna(ic2)*c(i2)+a_dna(ic3)*c(i3)

          mag = s*( alpha(i)+gamma(i)-phi_o(coreInd) )
          torque(i1) = torque(i1) + mag*ada
          torque(i2) = torque(i2) + mag*adb
          torque(i3) = torque(i3) + mag*adc
          mag = -s*( alpha(im1)+gamma(im1)-phi_o(coreInd) )
          torque(i1) = torque(i1) + mag*co
          torque(i2) = torque(i2) + mag*si
          torque(i3) = torque(i3) + 0.d0

c  extra Bending torques:
          torque(i2) = torque(i2) - g*beta_p(ic)*adc / dsin(beta_p(ic))
          torque(i3) = torque(i3) + g*beta_p(ic)*adb / dsin(beta_p(ic))

          ada = a(ib1)*a(i1)+a(ib2)*a(i2)+a(ib3)*a(i3)
          adb = a(ib1)*b(i1)+a(ib2)*b(i2)+a(ib3)*b(i3)
          adc = a(ib1)*c(i1)+a(ib2)*c(i2)+a(ib3)*c(i3)
          torque(i1) = torque(i1) + g*beta(im1)*
     +                 ( si*adc ) / dsin(beta(im1))
          torque(i2) = torque(i2) - g*beta(im1)*
     +                 ( co*adc ) / dsin(beta(im1))
          torque(i3) = torque(i3) + g*beta(im1)*
     +                 ( co*adb - si*ada ) / dsin(beta(im1))
c  extra Twisting torques:

          s1 = dsin( alpha_p(ic) )
          c1 = dcos( alpha_p(ic) )
          mag = s*(alpha(i)+gamma(i)-phi_o(coreInd))
     +          *dtan( 0.5d0*beta_p(ic) )
          cda = c_dna(ic1)*a(i1)+c_dna(ic2)*a(i2)+c_dna(ic3)*a(i3)
          bda = b_dna(ic1)*a(i1)+b_dna(ic2)*a(i2)+b_dna(ic3)*a(i3)
          cdb = c_dna(ic1)*b(i1)+c_dna(ic2)*b(i2)+c_dna(ic3)*b(i3)
          bdb = b_dna(ic1)*b(i1)+b_dna(ic2)*b(i2)+b_dna(ic3)*b(i3)
          cdc = c_dna(ic1)*c(i1)+c_dna(ic2)*c(i2)+c_dna(ic3)*c(i3)
          bdc = b_dna(ic1)*c(i1)+b_dna(ic2)*c(i2)+b_dna(ic3)*c(i3)
          torque(i1) = torque(i1) - mag*(s1*cda + c1*bda)
          torque(i2) = torque(i2) - mag*(s1*cdb + c1*bdb)
          torque(i3) = torque(i3) - mag*(s1*cdc + c1*bdc)

          s1 = dsin( gamma(im1) )
          c1 = dcos( gamma(im1) )
          mag = s*(alpha(im1)+gamma(im1)-phi_o(coreInd))*
     +          dtan( 0.5d0*beta(im1))
          cda = 0.d0
          bda = -si
          cdb = 0.d0
          bdb = co
          cdc = 1.0d0
          bdc = 0.d0
          torque(i1) = torque(i1) - mag*(s1*cda - c1*bda)
          torque(i2) = torque(i2) - mag*(s1*cdb - c1*bdb)
          torque(i3) = torque(i3) - mag*(s1*cdc - c1*bdc)

        end if
  200 continue
c
      i = n
      im1 = i-1
      i1 = 3*im1 + 1
      i2 = i1 + 1
      i3 = i2 + 1
      torque(i1) = -s*( alpha(im1)+gamma(im1)-phi_o(coreInd) )
      torque(i2) = 0.d0
      torque(i3) = 0.d0

!!!!!!!!!!!!!!!!!!!!!!!!! Display forces on file
      write(unit=199, fmt=1090) 'DNA Forces: '
      write(unit=199, fmt=1090) '1-Strecthing Forces: '

      write(unit=199, fmt=1090) 'f_DNA_stretch_x= [ '
      do 130 i = 1,n-1
          write(unit=199, fmt=1091) force_dnaStr_X(i) , ' , '
130   continue
          write(unit=199, fmt=1091) force_dnaStr_X(n) , ' ] '      
      
      write(unit=199, fmt=1090) 'f_DNA_stretch_y= [ '
      do 131 i = 1,n-1
          write(unit=199, fmt=1091) force_dnaStr_Y(i) , ' , '
131   continue
          write(unit=199, fmt=1091) force_dnaStr_Y(n) , ' ] '      
      
          write(unit=199, fmt=1090) 'f_DNA_stretch_z= [ ' 
      do 132 i = 1,n-1
          write(unit=199, fmt=1091) force_dnaStr_Z(i) , ' , '
132   continue
          write(unit=199, fmt=1091) force_dnaStr_Z(n) , ' ] '      
   
      
      write(unit=199, fmt=1090) '1-Bending Forces:'
      
      write(unit=199, fmt=1090) 'f_DNA_Bend_x= [ '
      do 133 i = 1,n-1
          write(unit=199, fmt=1091) force_dnaBend_X(i) , ' , '
133   continue
          write(unit=199, fmt=1091) force_dnaBend_X(n) , ' ] '      
      
      write(unit=199, fmt=1090) 'f_DNA_Bend_y= [ '
      do 134 i = 1,n-1
          write(unit=199, fmt=1091) force_dnaBend_Y(i) , ' , '
134   continue
          write(unit=199, fmt=1091) force_dnaBend_Y(n) , ' ] '      
      
          write(unit=199, fmt=1090) 'f_DNA_Bend_z= [ '
      do 135 i = 1,n-1
          write(unit=199, fmt=1091) force_dnaBend_Z(i) , ' , '
135   continue
          write(unit=199, fmt=1091) force_dnaBend_Z(n) , ' ] '      
   
   
      
      write(unit=199, fmt=1090) '1-Twisting Forces:'

         write(unit=199, fmt=1090) 'f_DNA_Twisting_x= [ '
      do 136 i = 1,n-1
          write(unit=199, fmt=1091) force_dnaTwist_X(i) , ' , '
136   continue
          write(unit=199, fmt=1091) force_dnaTwist_X(n) , ' ] '      
      
      write(unit=199, fmt=1090) 'f_DNA_Twisting_y= [ '
      do 137 i = 1,n-1
          write(unit=199, fmt=1091) force_dnaTwist_Y(i) , ' , '
137   continue
          write(unit=199, fmt=1091) force_dnaTwist_Y(n) , ' ] '      
      
          write(unit=199, fmt=1090) 'f_DNA_Twisting_z= [ '
      do 138 i = 1,n-1
          write(unit=199, fmt=1091) force_dnaTwist_Z(i) , ' , '
138   continue
          write(unit=199, fmt=1091) force_dnaTwist_Z(n) , ' ] '      

      
cccccccc       2nd     Part      - Abotaleb -                cccccccc
ccccccc      Mechanical forces on the tail beads cccccccc
ccccccc    internal [stretching and bending ] accounts only for the tail beads
ccccccc    external [stretching of the fixed beads with the parent core
ccccccc    accounts for both core and tail beads

C Divide up processors evenly 
          t_interval = t_n / np 
          t_par_start = myid * t_interval + 1
          t_par_stop = t_par_start + t_interval -1
          
      
      t_Eb = 0.0d0
      t_Ea = 0.0d0
c ... t_Eb = 1/2 * kb * (bi - b0)^2
c ... bond stretching force

c Stretching force on the tail beads [internal]

c
      do i = t_par_start, t_par_stop - 1
          
        if ( t_grp(i+1) .eq. t_grp(i) ) then
              
              x10 = t_X(i+1) - t_X(i)  
              y10 = t_Y(i+1) - t_Y(i) 
              z10 = t_Z(i+1) - t_Z(i)
              
            if(i .ne. 1) then
                if ( t_grp(i) .eq. t_grp(i-1) ) then
c
c (1.) Bending business:
c
                    Aim1(1) = Ai(1)
                    Aim1(2) = Ai(2)
                    Aim1(3) = Ai(3)
                    Bim1(1) = Bi(1)
                    Bim1(2) = Bi(2)
                    Bim1(3) = Bi(3)

                 
                x12 = t_X(i) - t_X(i-1)
                y12 = t_Y(i) - t_Y(i-1)
                z12 = t_Z(i) - t_Z(i-1)
          ! x12,y12,z12 now represent vector connecting tail bead #i with the prevoius one
                    if ( t_angle(i-1) .ge. 1.0d-10 ) then
                        g2 = (t_angle(i-1)-t_angle_v(i-1)) / 
     +                       ( dsin( t_angle(i-1) )*t_bond(i) )
                    else
                        g2 = 1.0d0 / t_bond(i)
                    end if
                  c2 = dcos( t_angle(i-1) )
                  Bi(1) = g2*( x12-c2*x10 )
                  Bi(2) = g2*( y12-c2*y10 )
                  Bi(3) = g2*( z12-c2*z10 )
                    
                    Strim1(1) = Stri(1)
                    Strim1(2) = Stri(2)
                    Strim1(3) = Stri(3)

                else
                    Aim1(1) = 0.d0
                    Aim1(2) = 0.d0
                    Aim1(3) = 0.d0
                
                    Bim1(1) = 0.d0
                    Bim1(2) = 0.d0
                    Bim1(3) = 0.d0

                    Bi(1) = 0.d0  ! if not precedded by a tail of same group
                    Bi(2) = 0.d0  ! Bi=0 [to be checked as in DNA linker/core model
                    Bi(3) = 0.d0  ! any linker DNA must be precedded by a core 
                  
                    Strim1(1) = 0.d0
                    Strim1(2) = 0.d0
                    Strim1(3) = 0.d0
                endif
            else
                    Aim1(1) = 0.d0
                    Aim1(2) = 0.d0
                    Aim1(3) = 0.d0
                    Bim1(1) = 0.d0
                    Bim1(2) = 0.d0
                    Bim1(3) = 0.d0

                    Bi(1) = 0.d0  ! if not precedded by a tail of same group
                    Bi(2) = 0.d0  ! Bi=0 [to be checked as in DNA linker/core model
                    Bi(3) = 0.d0  ! any linker DNA must be precedded by a core 
                    
                    Strim1(1) = 0.d0
                    Strim1(2) = 0.d0
                    Strim1(3) = 0.d0
            endif
            
          im1=i-1
          i1=3*im1+1
          i2=i1+1
          i3=i2+1

        if(i .le. t_par_stop - 2) then
            if ( t_grp(i+2) .eq. t_grp(i+1) ) then
              
               x12 = t_X(i+2) - t_X(i+1)
               y12 = t_Y(i+2) - t_Y(i+1)
               z12 = t_Z(i+2) - t_Z(i+1)
          
              if ( t_angle(i) .ge. 1.0d-10 ) then
                 g1 =(t_angle(i)-t_angle_v(i))
     +             /(dsin( t_angle(i) )*t_bond(i))
              else 
                 g1 = 1.0d0 / t_bond(i)
             end if
            
              c1 = dcos( t_angle(i) )
              Ai(1) = g1*( x12-c1*x10 )
              Ai(2) = g1*( y12-c1*y10 )
              Ai(3) = g1*( z12-c1*z10 )
          else
            Ai(1) = 0.d0
            Ai(2) = 0.d0
            Ai(3) = 0.d0
          endif
        else
            Ai(1) = 0.d0
            Ai(2) = 0.d0
            Ai(3) = 0.d0
        endif
        
            
            
          Stri(1) = (t_bond(i)-t_bond_v(i))*x10  ! check that a direction is the direction connecting consecutive beads
          Stri(2) = (t_bond(i)-t_bond_v(i))*y10  ! update_mod recovers a,b,c get_euler doesn't  
          Stri(3) = (t_bond(i)-t_bond_v(i))*z10  ! a actually in the direction connecting two consecutive beads 

c forward segment...
        force_tStrInt_X(i) = t_bond_c(i)*( Stri(1) )
        force_tStrInt_Y(i) = t_bond_c(i)*( Stri(2) )
        force_tStrInt_Z(i) = t_bond_c(i)*( Stri(3) )
        
        force_tBend_X(i)= - t_angle_c(i)*( Ai(1)+Bi(1) )
        force_tBend_Y(i)= - t_angle_c(i)*( Ai(2)+Bi(2) )
        force_tBend_Z(i)= - t_angle_c(i)*( Ai(3)+Bi(3) )
        
        force_t(i1) = force_tStrInt_X(i) + force_tBend_X(i)
        force_t(i2) = force_tStrInt_Y(i) + force_tBend_Y(i)
        force_t(i3) = force_tStrInt_Z(i) + force_tBend_Z(i)
       
        
c back segment...
        force_tStrInt_X(i) =force_tStrInt_X(i)+ t_bond_c(i)*(-Strim1(1))
        force_tStrInt_Y(i) =force_tStrInt_Y(i)+ t_bond_c(i)*(-Strim1(2))
        force_tStrInt_Z(i) =force_tStrInt_Z(i)+ t_bond_c(i)*(-Strim1(3))
        
        force_tBend_X(i)=force_tBend_X(i)-t_angle_c(i)
     +                   *(-Aim1(1)-Bim1(1) )
        force_tBend_Y(i)=force_tBend_Y(i)-t_angle_c(i)
     +                   *(-Aim1(2)-Bim1(2) )
        force_tBend_Z(i)=force_tBend_Z(i)-t_angle_c(i)
     +                   *(-Aim1(3)-Bim1(3) )
        
        force_t(i1) =force_t(i1) + force_tStrInt_X(i) + force_tBend_X(i)
        force_t(i2) =force_t(i2) + force_tStrInt_Y(i) + force_tBend_Y(i)
        force_t(i3) =force_t(i3) + force_tStrInt_Z(i) + force_tBend_Z(i)
        
        endif
                end do
                
                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
c  Numerical Integration to get the energy [checking]
c  t_Eb_num,t_Ea_num

        ! Now it is done using trapezoidal unqual spaced
        ! for a better performance search for simpson unequal spaced

ccccccccccccccccc  strecthing 

 !          call Isimpson(t_n,t_X,force_tStrInt_X,temp)
 !          t_Eb_num=temp;
 !          call Isimpson(t_n,t_Y,force_tStrInt_Y,temp)
 !          t_Eb_num=t_Eb_num+temp;
 !          call Isimpson(t_n,t_Z,force_tStrInt_Z,temp)
 !          t_Eb_num=t_Eb_num+temp;
 !          t_Eb_num=-t_Eb_num; ! E = - Integral F .dl  in calculating the derivative original work dropped the minus sign
cccccccccccccccc  bending
        
 !          call Isimpson(t_n,t_X,force_tBend_X,temp)
 !          t_Ea_num=temp;
 !          call Isimpson(t_n,t_Y,force_tBend_Y,temp)
 !          t_Ea_num=t_Eb_num+temp;
 !          call Isimpson(t_n,t_Z,force_tBend_Z,temp)
 !          t_Ea_num=t_Eb_num+temp;
 !          t_Ea_num=t_Ea_num! E = - Integral F .dl  
              

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Electrostatic parameters:
      ql_ql = q_l * q_l
c Exclude volume parameters:
c
c (2) Bead-Bead interactions: (a.) electrostatic calculations
c
      do 300 i = 1,nm1
        i1 = 3*(i-1) + 1
        i2 = i1 + 1
        i3 = i2 + 1

        do 300 j = (i+1),n
          j1 = 3*(j-1) + 1
          j2 = j1 + 1
          j3 = j2 + 1
          dist = dsqrt( ( r(j1)-r(i1) )**2 +
     +                  ( r(j2)-r(i2) )**2 +
     +                  ( r(j3)-r(i3) )**2 )

          if ( dist .le. Rcut ) then
          if ( type(i) .eq. 0 ) then
            if ( type(j). eq. 0 ) then
              if ( abs(i-j) .gt. 1 ) then
c
c Linker-Linker interaction:
c
c                dist = dsqrt( ( r(j1)-r(i1) )**2 + 
c     +                        ( r(j2)-r(i2) )**2 + 
c     +                        ( r(j3)-r(i3) )**2 )
                mi = 1.0d0/dist
                z(1) = mi*( r(i1) - r(j1) )
                z(2) = mi*( r(i2) - r(j2) )
                z(3) = mi*( r(i3) - r(j3) )
                g1 = k_e*ql_ql*dexp(-debye*dist)*(debye*dist+1.d0)
     +                    / (dist*dist)
                force(i1) = force(i1) + g1*z(1)
                force(i2) = force(i2) + g1*z(2)
                force(i3) = force(i3) + g1*z(3)
                force(j1) = force(j1) - g1*z(1)
                force(j2) = force(j2) - g1*z(2)
                force(j3) = force(j3) - g1*z(3)
              end if

            else
c
c Linker-Core interaction ( i is a linker; j is a core ):
c 
            do k = 1,Nq
              k1 = 3*(k-1) + 1
              k2 = k1 + 1
              k3 = k2 + 1
              z(1) = ( r(i1) - (r(j1)+a(j1)*core_pos(k1)
     +                               +b(j1)*core_pos(k2)
     +                               +c(j1)*core_pos(k3)) )
              z(2) = ( r(i2) - (r(j2)+a(j2)*core_pos(k1)
     +                               +b(j2)*core_pos(k2)
     +                               +c(j2)*core_pos(k3)) )
              z(3) = ( r(i3) - (r(j3)+a(j3)*core_pos(k1)
     +                               +b(j3)*core_pos(k2)
     +                               +c(j3)*core_pos(k3)) )
              dist = dsqrt( z(1)**2 + z(2)**2 + z(3)**2 )
              mi = 1.0d0/dist
              z(1) = mi*z(1)
              z(2) = mi*z(2)
              z(3) = mi*z(3)
              if ( abs(i-j) .gt. 1 ) then
                g1 = k_e*q_l*core_q(k)
     +               *dexp(-debye*dist)*(debye*dist+1.d0)
     +                    / (dist*dist)
                force(i1) = force(i1) + g1*z(1)
                force(i2) = force(i2) + g1*z(2)
                force(i3) = force(i3) + g1*z(3)
                force(j1) = force(j1) - g1*z(1)
                force(j2) = force(j2) - g1*z(2)
                force(j3) = force(j3) - g1*z(3)
c  torque due to e forces on j-core :
                fa = -g1*( a(j1)*z(1)+a(j2)*z(2)+a(j3)*z(3) )
                fb = -g1*( b(j1)*z(1)+b(j2)*z(2)+b(j3)*z(3) ) 
                fc = -g1*( c(j1)*z(1)+c(j2)*z(2)+c(j3)*z(3) )
                torque(j1) = torque(j1) +
     +            fc*core_pos(k2) - fb*core_pos(k3)
                torque(j2) = torque(j2) +
     +            fa*core_pos(k3) - fc*core_pos(k1)
                torque(j3) = torque(j3) +
     +            fb*core_pos(k1) - fa*core_pos(k2)
              end if
c
c Excluded volume forces:
c
              if ((dist .le. 8.0d0) .and. (core_q(k) .gt. 0.d0)) then
                s1 = 3.0d0
                s2 = 3.0d0
                g1 = k_ex*( (12.0d0/s1)*(s1/dist)**13
     +                     -( 6.0d0/s2)*(s2/dist)**7)
                force(i1) = force(i1) + g1*z(1)
                force(i2) = force(i2) + g1*z(2)
                force(i3) = force(i3) + g1*z(3)
                force(j1) = force(j1) - g1*z(1)
                force(j2) = force(j2) - g1*z(2)
                force(j3) = force(j3) - g1*z(3)
                fa = -g1*( a(j1)*z(1)+a(j2)*z(2)+a(j3)*z(3) )
                fb = -g1*( b(j1)*z(1)+b(j2)*z(2)+b(j3)*z(3) )
                fc = -g1*( c(j1)*z(1)+c(j2)*z(2)+c(j3)*z(3) )
                torque(j1) = torque(j1) +
     +            fc*core_pos(k2) - fb*core_pos(k3)
                torque(j2) = torque(j2) +
     +            fa*core_pos(k3) - fc*core_pos(k1)
                torque(j3) = torque(j3) +
     +            fb*core_pos(k1) - fa*core_pos(k2)
              end if

            end do

            end if

          else
            if ( type(j). eq. 0 ) then
c
c Core-Linker interaction ( i is a core; j is a linker ):
c
              do k = 1,Nq
                k1 = 3*(k-1) + 1
                k2 = k1 + 1
                k3 = k2 + 1
                z(1) = ( -r(j1) + (r(i1)+a(i1)*core_pos(k1)
     +                                  +b(i1)*core_pos(k2)
     +                                  +c(i1)*core_pos(k3)) )
                z(2) = ( -r(j2) + (r(i2)+a(i2)*core_pos(k1)
     +                                  +b(i2)*core_pos(k2)
     +                                  +c(i2)*core_pos(k3)) )
                z(3) = ( -r(j3) + (r(i3)+a(i3)*core_pos(k1)
     +                                  +b(i3)*core_pos(k2)
     +                                  +c(i3)*core_pos(k3)) )
		dist = dsqrt( z(1)**2 + z(2)**2 + z(3)**2 )

		mi = 1.0d0/dist
		z(1) = mi*z(1)
		z(2) = mi*z(2)
		z(3) = mi*z(3)
                if ( abs(i-j) .gt. 1 ) then
		  g1 = k_e*q_l*core_q(k)
     +                 *dexp(-debye*dist)*(debye*dist+1.d0)
     +                 / (dist*dist)
                  force(i1) = force(i1) + g1*z(1)
                  force(i2) = force(i2) + g1*z(2)
                  force(i3) = force(i3) + g1*z(3)
                  force(j1) = force(j1) - g1*z(1)
                  force(j2) = force(j2) - g1*z(2)
                  force(j3) = force(j3) - g1*z(3)
c  torque due to e forces on i-core :
                  fa = g1*( a(i1)*z(1)+a(i2)*z(2)+a(i3)*z(3) )
                  fb = g1*( b(i1)*z(1)+b(i2)*z(2)+b(i3)*z(3) ) 
                  fc = g1*( c(i1)*z(1)+c(i2)*z(2)+c(i3)*z(3) )
                  torque(i1) = torque(i1) +
     +              fc*core_pos(k2) - fb*core_pos(k3) 
                  torque(i2) = torque(i2) +
     +              fa*core_pos(k3) - fc*core_pos(k1) 
                  torque(i3) = torque(i3) +
     +              fb*core_pos(k1) - fa*core_pos(k2)
                end if
c
c Excluded volume forces:
c
                if ((dist .le. 8.0d0) .and. (core_q(k) .gt. 0.d0)) then
                  s1 = 3.0d0
                  s2 = 3.0d0
                  g1 = k_ex*( (12.0d0/s1)*(s1/dist)**13
     +                       -( 6.0d0/s2)*(s2/dist)**7)
                  force(i1) = force(i1) + g1*z(1)
                  force(i2) = force(i2) + g1*z(2)
                  force(i3) = force(i3) + g1*z(3)
                  force(j1) = force(j1) - g1*z(1)
                  force(j2) = force(j2) - g1*z(2)
                  force(j3) = force(j3) - g1*z(3)
                  fa = g1*( a(i1)*z(1)+a(i2)*z(2)+a(i3)*z(3) )
                  fb = g1*( b(i1)*z(1)+b(i2)*z(2)+b(i3)*z(3) )
                  fc = g1*( c(i1)*z(1)+c(i2)*z(2)+c(i3)*z(3) )
                  torque(i1) = torque(i1) +
     +              fc*core_pos(k2) - fb*core_pos(k3)
                  torque(i2) = torque(i2) +
     +              fa*core_pos(k3) - fc*core_pos(k1)
                  torque(i3) = torque(i3) +
     +              fb*core_pos(k1) - fa*core_pos(k2)
                end if

	      end do

	    else
c
c Core-Core interaction:
c
	      do k = 1,Nq
		k1 = 3*(k-1)+1
		k2 = k1+1
		k3 = k2+1
		do l = 1,Nq
		  l1 = 3*(l-1)+1
		  l2 = l1+1
		  l3 = l2+1
                  z(1) = r(i1)+a(i1)*core_pos(k1) +
     +                         b(i1)*core_pos(k2) +
     +                         c(i1)*core_pos(k3) -
     +                   r(j1)-a(j1)*core_pos(l1) -
     +                         b(j1)*core_pos(l2) -
     +                         c(j1)*core_pos(l3)
                  z(2) = r(i2)+a(i2)*core_pos(k1) +
     +                         b(i2)*core_pos(k2) +
     +                         c(i2)*core_pos(k3) -
     +                   r(j2)-a(j2)*core_pos(l1) -
     +                         b(j2)*core_pos(l2) -
     +                         c(j2)*core_pos(l3) 
                  z(3) = r(i3)+a(i3)*core_pos(k1) +
     +                         b(i3)*core_pos(k2) +
     +                         c(i3)*core_pos(k3) -
     +                   r(j3)-a(j3)*core_pos(l1) -
     +                         b(j3)*core_pos(l2) -
     +                         c(j3)*core_pos(l3) 

                  dist = dsqrt( z(1)**2 + z(2)**2 + z(3)**2 )
                  mi = 1.0d0/dist
                  z(1) = mi*z(1)
                  z(2) = mi*z(2)
                  z(3) = mi*z(3)
                  g1 = k_e*core_q(k)*core_q(l)
     +                 *dexp(-debye*dist)*(debye*dist+1.d0)
     +                    / (dist*dist)
                  force(i1) = force(i1) + g1*z(1)
                  force(i2) = force(i2) + g1*z(2)
                  force(i3) = force(i3) + g1*z(3)
                  force(j1) = force(j1) - g1*z(1)
                  force(j2) = force(j2) - g1*z(2)
                  force(j3) = force(j3) - g1*z(3)
c  torque due to e forces on i-core :
                  fa = g1*( a(i1)*z(1)+a(i2)*z(2)+a(i3)*z(3) )
                  fb = g1*( b(i1)*z(1)+b(i2)*z(2)+b(i3)*z(3) ) 
                  fc = g1*( c(i1)*z(1)+c(i2)*z(2)+c(i3)*z(3) )
                  torque(i1) = torque(i1) +
     +              fc*core_pos(k2) - fb*core_pos(k3)
                  torque(i2) = torque(i2) +
     +              fa*core_pos(k3) - fc*core_pos(k1)
                  torque(i3) = torque(i3) +
     +              fb*core_pos(k1) - fa*core_pos(k2)
c  torque due to e forces on j-core :
                  fa = -g1*( a(j1)*z(1)+a(j2)*z(2)+a(j3)*z(3) )
                  fb = -g1*( b(j1)*z(1)+b(j2)*z(2)+b(j3)*z(3) )
                  fc = -g1*( c(j1)*z(1)+c(j2)*z(2)+c(j3)*z(3) )
                  torque(j1) = torque(j1) +
     +              fc*core_pos(l2) - fb*core_pos(l3)
                  torque(j2) = torque(j2) +
     +              fa*core_pos(l3) - fc*core_pos(l1)
                  torque(j3) = torque(j3) +
     +              fb*core_pos(l1) - fa*core_pos(l2)
c
c Excluded volume forces:
c
                  if ( dist .le. 8.0d0 ) then
                    s1 = 2.0d0
                    s2 = 2.0d0
                    g1 = k_ex*( (12.0d0/s1)*(s1/dist)**13
     +                         -( 6.0d0/s2)*(s2/dist)**7)
                    force(i1) = force(i1) + g1*z(1)
                    force(i2) = force(i2) + g1*z(2)
                    force(i3) = force(i3) + g1*z(3)
                    force(j1) = force(j1) - g1*z(1)
                    force(j2) = force(j2) - g1*z(2)
                    force(j3) = force(j3) - g1*z(3)
                    fa = g1*( a(i1)*z(1)+a(i2)*z(2)+a(i3)*z(3) )
                    fb = g1*( b(i1)*z(1)+b(i2)*z(2)+b(i3)*z(3) )
                    fc = g1*( c(i1)*z(1)+c(i2)*z(2)+c(i3)*z(3) )
                    torque(i1) = torque(i1) +
     +                fc*core_pos(k2) - fb*core_pos(k3)
                    torque(i2) = torque(i2) +
     +                fa*core_pos(k3) - fc*core_pos(k1)
                    torque(i3) = torque(i3) +
     +                fb*core_pos(k1) - fa*core_pos(k2)
                    fa = -g1*( a(j1)*z(1)+a(j2)*z(2)+a(j3)*z(3) )
                    fb = -g1*( b(j1)*z(1)+b(j2)*z(2)+b(j3)*z(3) )
                    fc = -g1*( c(j1)*z(1)+c(j2)*z(2)+c(j3)*z(3) )
                    torque(j1) = torque(j1) +
     +                fc*core_pos(l2) - fb*core_pos(l3)
                    torque(j2) = torque(j2) +
     +                fa*core_pos(l3) - fc*core_pos(l1)
                    torque(j3) = torque(j3) +
     +                fb*core_pos(l1) - fa*core_pos(l2)
                  end if

              end do

            end do
            end if

          end if
c end of Bead-Bead cutoff
        end if 

  300 continue
       close(unit=199)
!!!!!!!!!!!!!!!!!  testing
!!!!! testing the derivative using testgh routine
! test the stretching energy derivative
      xc=r
      yhy=0
      num=n
      isAnalytic=1
        aold=a
        bold=b
        cold=c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! testing strecthing energy 
            fnc=E(1)
       do 1334 i = 1,n

          i1 = 3*(i-1) + 1
          i2 = i1 + 1
          i3 = i2 + 1
          y(i1)=0
          y(i2)=0
          y(i3)=0
 
          gc(i1)=-force_dnaStr_X(i)
          gc(i2)=-force_dnaStr_Y(i)
          gc(i3)=-force_dnaStr_Z(i)
1334    continue
       y(1)=0.05
       
      call       testgh(n3,xc,fnc,gc,y,yhy,vec,
     +        n_c,nc3,num,n3,type,r,ro,d1,go,a,b,c,
     +        aold,bold,cold,
     +        alpha,beta,gamma,length,a_dna,b_dna,c_dna,
     +        alpha_p,beta_p,gamma_p,
     +    lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,E,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    1,isAnalytic,
     +    np,myid,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! testing bending energy 
            fnc=E(2)
       do 1335 i = 1,n

          i1 = 3*(i-1) + 1
          i2 = i1 + 1
          i3 = i2 + 1
          y(i1)=0
          y(i2)=0
          y(i3)=0
      
          gc(i1)=-force_dnaBend_X(i)
          gc(i2)=-force_dnaBend_Y(i)
          gc(i3)=-force_dnaBend_Z(i)
             
          
1335    continue
       y(1)=0.05
      call       testgh(n3,xc,fnc,gc,y,yhy,vec,
     +        n_c,nc3,num,n3,type,r,ro,d1,go,a,b,c,
     +        aold,bold,cold,
     +        alpha,beta,gamma,length,a_dna,b_dna,c_dna,
     +        alpha_p,beta_p,gamma_p,
     +    lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,E,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    2,isAnalytic,
     +    np,myid,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! testing twisting energy 
            fnc=E(3)
       do 1336 i = 1,n

          i1 = 3*(i-1) + 1
          i2 = i1 + 1
          i3 = i2 + 1
          y(i1)=0
          y(i2)=0
          y(i3)=0

          gc(i1)=-force_dnaTwist_X(i)
          gc(i2)=-force_dnaTwist_Y(i)
          gc(i3)=-force_dnaTwist_Z(i)
      
          
1336    continue
       y(1)=0.05
      call       testgh(n3,xc,fnc,gc,y,yhy,vec,
     +        n_c,nc3,num,n3,type,r,ro,d1,go,a,b,c,
     +        aold,bold,cold,
     +        alpha,beta,gamma,length,a_dna,b_dna,c_dna,
     +        alpha_p,beta_p,gamma_p,
     +    lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,E,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    3,isAnalytic,
     +    np,myid,ierr)

 1090 format(1x,A45)

 1091 format(1x,F12.9, 2x,A5)

      end



