cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine potential(n_c,n,n3,type,r,a,b,c,alpha,beta,gamma,
     +	  length,beta_p,lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,E,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    np,myid,ierr,fcount)

      use modglob
      use mpi          
      implicit NONE

      integer n_c,n,n3, type(n), h_n
      double precision r(n3), a(n3),b(n3),c(n3)
      double precision alpha(n),beta(n),gamma(n), length(n)
      double precision beta_p(n_c)
      double precision lo,hd2,gd2,sd2,k_e,k_ex
      integer Nq, Nq3,e_counter
      double precision core_pos(Nq3), core_q(Nq)
      double precision debye, q_l, phi_o(n_c)
      double precision E(6),Esum(6)
      double precision pot_send_array(28),sum_pot_send_array(28)

      integer i, i1,i2,i3, j, j1,j2,j3, ib
      integer k, k1,k2,k3, l, l1,l2,l3
      integer nm1
      double precision dist, z(3)
      double precision Ev,Ec,ql_ql, ssle,ssbe,sstw
      double precision s1, s2, Rcut, Rcut2

      integer t_n
      double precision t_chg(t_n), h_chg(h_n)
      double precision t_X(t_n), t_Y(t_n), t_Z(t_n)

      double precision h_X(h_n), h_Y(h_n), h_Z(h_n)
      double precision t_X0(t_n), t_Y0(t_n), t_Z0(t_n)
      double precision t_bond(t_n), t_bond_v(t_n), t_bond_c(t_n)
      double precision t_angle(t_n), t_angle_v(t_n), t_angle_c(t_n)
      double precision t_exv_e,t_exv_d,evd_hcG,evd_hcC,evd_hlG,evd_hlC
      double precision t_Eb, t_Ea, t_Ec, t_Ev
      integer t_grp(t_n), t_fix(t_n)
      double precision vdw_cut, vdw_cut2
      double precision evd_tc, evd_tl, evd_cc, evd_cl, evd_ll, evd_hh
      integer m1, m2, m3
      double precision dx, dy, dz, dr2, dr, dt2
      double precision x12, y12, z12, r12, r12s
      double precision x10, y10, z10, r10, r10s
      double precision c0, phi
      double precision test1, test2, test3
      double precision SMALL, SMALL2
      parameter(SMALL = 1.0d-4, SMALL2 = 1.0d-8)
c ... The 6 components of electrostatics, L: linker; C: core; T: tail
      double precision Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL, h_tc
c ... The 2 subcomponents of Ec_TT: Ec_TT1: tails in the same core; Ec_TT2: tails in different cores
c     The 2 subcomponents of Ec_CT: Ec_CT1: tails with its own core; Ec_CT2: tails with other cores 
      double precision Ec_TT1, Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH
      integer nt, ntpc
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
      
      integer index,indexj(n)

c     parallelization variables   
        integer np,myid,ierr,n_interval,t_interval,h_interval
        integer n_par_start,n_par_stop,t_par_start,t_par_stop
        integer h_par_start,h_par_stop  
        integer fcount
cccccccccccccccccccccccccccccccccccccccccccc
  ! Abotaleb  16 Nov      
cccccccccccccccccccccccccccccccccccccccccccc

        character(len=1024) :: filename
      write (filename, "(A20,I3,A5)") "Energies data ", fcount,'.txt'

          open(unit=133,name=filename,access='SEQUENTIAL',
     +        status='unknown')
        
      index=0
      do i=1,n
         if(type(i).ne.0) index=index+1
         indexj(i)=index
      enddo



c     number of beads (core and DNA) in each repeating unit (monomer)
      nt  = n / n_c
c     number of tail beads per core.
      ntpc = t_n/n_c
      nm1 = n-1
      Rcut2 = Rcut * Rcut
      vdw_cut2 = vdw_cut * vdw_cut
    

C Divide up processors evenly 
          t_interval = t_n / np 
          t_par_start = myid * t_interval + 1
          t_par_stop = t_par_start + t_interval -1
          
          n_interval = n / np
          n_par_start = myid *  n_interval +1
          n_par_stop = n_par_start + n_interval -1

          h_interval = h_n / np
          h_par_start = myid *  h_interval +1
          h_par_stop = h_par_start + h_interval -1
          
          if (myid .eq. np-1) then
          t_par_stop = t_n
          n_par_stop = n
          h_par_stop = h_n
          endif


c     initialize
      Ec = 0.d0
      Ev = 0.d0
      Ec_LL = 0.0d0
      Ec_CC = 0.0d0
      Ec_TT = 0.0d0  
      Ec_TT1 = 0.0d0
      Ec_TT2 = 0.0d0
      Ec_LC = 0.0d0
      Ec_CT = 0.0d0
      Ec_CT1 = 0.0d0
      Ec_CT2 = 0.0d0
      Ec_TL = 0.0d0
      t_Ec = 0.0d0
      t_Ev = 0.0d0
      Eb_tc = 0.0d0
      Ev_LL = 0.0d0
      Ev_CC = 0.0d0
      Ev_TT = 0.0d0
      Ev_LC = 0.0d0
      Ev_CT = 0.0d0
      Ev_TL = 0.0d0
      Ev_TT1 = 0.0d0
      Ev_TT2 = 0.0d0
      Ev_CT1 = 0.0d0
      Ev_CT2 = 0.0d0
      Ev_HH = 0.0d0
      Ev_HT = 0.0d0
      Ev_HL1 = 0.0d0
      Ev_HL2 = 0.0d0
      Ev_HC = 0.0d0
cccccccccccccccccccccc
C	print*,'Starting mechanical potential'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Part 1: mechanical potential
c - Es/ssle: stretching
c - Eb/ssbe: bending
c - Et/sstw: twisting
      ssle = 0.d0
      ssbe = 0.d0
      sstw = 0.d0
C      if (myid.eq.0) then  
      write(unit=133, fmt=1089) 'DNA Twist According to Linker length:'
      do 100 i = 1,nm1
        ssle = ssle + ( length(i) - lo )**2
        ssbe = ssbe + ( beta(i) )**2
        ib=indexj(i)
        sstw = sstw + ( alpha(i)+gamma(i)-phi_o(ib) )**2
                
        write(unit=133, fmt=2089) 'for bead #',i,' ; index = ',ib       
        write(unit=133, fmt=2289) 'phi_o(',ib,') = ',phi_o(ib)
  100 continue
      
      do i = 1,n_c
        ssbe = ssbe + ( beta_p(i) )**2
      end do
      
      ssle = hd2*ssle 
      ssbe = gd2*ssbe
      sstw = sd2*sstw
C      endif  
c ********  Abotaleb
      write(unit=133, fmt=1089)'for DNA beads:'
              
          write(unit=133, fmt=1809)" 1-Stretching Energy = ",ssle
          write(unit=133, fmt=1809)"  -Stretching Const = " ,hd2
          write(unit=133, fmt=1809)" 2-Bending Energy = ",ssbe
          write(unit=133, fmt=1809)"  -Bending Const = " ,gd2
          write(unit=133, fmt=1809)" 3-Twisting Energy = ",sstw
          write(unit=133, fmt=1809)"  -Twisting Const = " ,sd2
          
          

C	print*,'Starting bead-bead pot'

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c Part 2: Bead-Bead (LL, LC, CL, CC) interactions (L: linker, C: core)
c - Ec: electrostatic calculations
c   Ec_LL  Ec_LC Ec_CC
c - Ev: excluded volume calcs.

      ql_ql = q_l * q_l

      do 200 i = n_par_start,n_par_stop-1
C     do 200 i = 1,nm1    
        i1 = 3*(i-1) + 1
        i2 = i1 + 1
        i3 = i2 + 1

        do 200 j = (i+1),n
          j1 = 3*(j-1) + 1
          j2 = j1 + 1
          j3 = j2 + 1
          dist = dsqrt( ( r(j1)-r(i1) )**2 +
     +                  ( r(j2)-r(i2) )**2 +
     +                  ( r(j3)-r(i3) )**2 ) 

          if (dist .le. Rcut+11.0) then  !! RCUT + approx 2*radius of core

          if ( type(i) .eq. 0 ) then
             if ( type(j). eq. 0 ) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ## 1. LL: Linker-Linker interaction:
c ..... LL/Ec: 
                if (dist .le. Rcut) then !!RCUT
                  if ( abs(i-j) .gt. 1 ) then  !<<<<<MADE IT 2 from 1
                     Ec_LL = Ec_LL + k_e*ql_ql*dexp(-debyell*dist)/dist                     
                   end if
        				
		        endif
           else
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ## 2. LC: Linker-Core interaction ( i is a linker; j is a core ):
               do k = 1,Nq
                  k1 = 3*(k-1) + 1
                  k2 = k1 + 1
                  k3 = k2 + 1
                  z(1) = r(i1) - ( r(j1)+a(j1)*core_pos(k1)
     +                              +b(j1)*core_pos(k2)
     +                              +c(j1)*core_pos(k3) ) 
                  z(2) = r(i2) - ( r(j2)+a(j2)*core_pos(k1)
     +                              +b(j2)*core_pos(k2)
     +                              +c(j2)*core_pos(k3) )
                  z(3) = r(i3) - ( r(j3)+a(j3)*core_pos(k1)
     +                              +b(j3)*core_pos(k2)
     +                              +c(j3)*core_pos(k3) )
                  dist = dsqrt( z(1)**2 + z(2)**2 + z(3)**2 )
c ..... LC/Ec:
                  if (dist .le. Rcut) then    
                    if ( abs(i-j) .gt. 1 ) then
                      Ec_LC = Ec_LC + k_e * 
     +                  q_l*core_q(k)*dexp(-debye*dist) / dist
                    end if
		          endif
c ..... LC/Ev:
                  if ( (abs(i-j).gt.1).and.(dist.le.vdw_cut)) then
                    s1 = evd_cl
                    s2 = evd_cl

                    Ev_LC = Ev_LC + k_ex * ((s1/dist)**12-(s2/dist)**6)

                  end if
               end do
	       
             end if
          else
            if ( type(j). eq. 0 ) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ## 3. CL: Core-Linker interaction ( i is a core; j is a linker ):
              do k = 1,Nq
                 k1 = 3*(k-1) + 1
                 k2 = k1 + 1
                 k3 = k2 + 1
                 z(1) = -r(j1) + ( r(i1)+a(i1)*core_pos(k1)
     +                                 +b(i1)*core_pos(k2)
     +                                 +c(i1)*core_pos(k3) )
                 z(2) = -r(j2) + ( r(i2)+a(i2)*core_pos(k1)
     +                                 +b(i2)*core_pos(k2)
     +                                 +c(i2)*core_pos(k3) )
                 z(3) = -r(j3) + ( r(i3)+a(i3)*core_pos(k1)
     +                                 +b(i3)*core_pos(k2)
     +                                 +c(i3)*core_pos(k3) )
                 dist = dsqrt( z(1)**2 + z(2)**2 + z(3)**2 )
c ..... CL/Ec:
                 if (dist .le. Rcut) then   
                 if ( abs(i-j) .gt. 1 ) then
                     Ec_LC = Ec_LC + k_e * 
     +                    q_l*core_q(k)*dexp(-debye*dist) / dist
                 end if
		 endif
c ..... CL/Ev:
                 if ( (abs(i-j) .gt. 1) .and. (dist .le. vdw_cut) ) then
c                if (dist .le. vdw_cut)  then
                   s1 = evd_cl
                   s2 = evd_cl

                    Ev_LC = Ev_LC + k_ex * ((s1/dist)**12-(s2/dist)**6)


                 end if
              end do
            else
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ## 4. Core-Core interaction:
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
c ..... CC/Ec:
                    if (dist .le. Rcut) then   
                       Ec_CC = Ec_CC + k_e * 
     +                 core_q(k)*core_q(l)*dexp(-debye*dist)/dist
                    endif
c ..... CC/Ev:
                    if ( dist .le. vdw_cut ) then
                      s1 = evd_cc
                      s2 = evd_cc

                      Ev_CC = Ev_CC + 
     +                     k_ex*((s1/dist)**12-(s2/dist)**6)

                      
                    end if
                end do
              end do

            end if
         end if
c End of cutoff-if
         end if  !! RCUT + core radii
  200 continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Part 2 Mechnical Energies Calculations
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C       print*,'starting LH mech pot' 
C       print*,'h_n',h_n
ccccccccccc
ccc Modification: TONI (LHref)
ccc linker histone - linker histone interactions (elastic force-field)
      ELHel = 0.
      ELHstr = 0.
      ELHben = 0.
      ELHtor = 0.
C      if(myid.eq.0) then  
C      if (withlink) then !................................WITHLINK IF START
      do i = h_par_start,h_par_stop
C     do i = 1,h_n  
         if(LHboundbds(i)) then
            !!! Stretching
            nconn = h_conn(i)
            if(nconn==1) then
               ktmp = rkstrd2(i) !!! Elastic constant
               dr = lengthLH(i) !!! Distance (next bead)
               dr = dr - rkstreq(i) !!! Relative to equilibrium distance
               dr2 = dr*dr
               e_tmp = ktmp*dr2 !!! Elastic term
               ELHstr = ELHstr + e_tmp
               !     write(*,*) 'ELHstr',i,ELHstr,lengthLH(i)        
            endif
            !!!   Bending
            nbend = h_bend(i)
            if(nbend==1) then
               ktmp = rkbend2(i) !!! Elastic constant
               dr = betaLH(i)   !!! Angle
               dr = dr - rkbeneq(i) !!! Relative to equilibrium angle
               dr2 = dr*dr
               e_tmp = ktmp*dr2 !!! Elastic term
               ELHben = ELHben + e_tmp
               !     write(*,*) 'ELHben',i,ELHben,dr            
            endif

         endif !!! LHboundbds(i)
      enddo
      ELHel = ELHstr + ELHben + ELHtor
      !write(*,*) 'ELHel',ELHel
C      endif  
c     linker histone - linker histone interactions
      Ec_HH = 0.0

      do i = h_par_start, h_par_stop-1
C     do i = 1,h_n-1   
        !!! Calculate the energy if the LH bead exists
         if(LHboundbds(i)) then
            !do i = 1, h_n
            !jlow = nbLH*((i-1)/nbLH+1)
            jlow = i - 1        !!! for same LH avoid near-neigh Ec Ev
            jhigh = i + 1
            grpi = LH_grp(i)
            do j = i+1, h_n
               !!! Calculate the energy if the LH bead exists
               if(LHboundbds(j)) then
cccccccccccccccccccccccc
                  grpj = LH_grp(j)
                  !!!   Different LHs: Ec and Ev are active
                  if(grpi.ne.grpj) then
                     dx = h_X(i) - h_X(j)
                     dy = h_Y(i) - h_Y(j)
                     dz = h_Z(i) - h_Z(j)
                     dr2 = dx*dx + dy*dy + dz*dz
                     if (dr2 .le. Rcut2) then
                        if (dr2 .gt. SMALL2) then
                           dr = sqrt(dr2)
                        else
                           dr = SMALL
                        end if
                        Ec_HH = Ec_HH + k_e*h_chg(i)*
     +                       h_chg(j)*dexp(-debye*dr)/dr
                     end if

                     dr = sqrt(dr2)

                     if(dr .le. vdw_cut) then
                        revd = (revd_hh(i) + revd_hh(j))/2
                        s1 = revd
                        s2 = revd
                        Ev_HH = Ev_HH + 
     +                       k_ex*((s1/dr)**12-(s2/dr)**6)
                     endif   
                  !!!   Same LH
                  !!!   Near-neighbors: Ec and Ev are inactive
                  !!!   Rest: Ev active   
                  else if((j.lt.jlow).OR.(j.gt.jhigh)) then
                     dx = h_X(i) - h_X(j)
                     dy = h_Y(i) - h_Y(j)
                     dz = h_Z(i) - h_Z(j)
                     dr2 = dx*dx + dy*dy + dz*dz
                     dr = sqrt(dr2)
                     if(dr .le. vdw_cut) then
                        revd = (revd_hh(i) + revd_hh(j))/2
                        s1 = revd
                        s2 = revd
                        Ev_HH = Ev_HH + 
     +                       k_ex*((s1/dr)**12-(s2/dr)**6)
                     endif   
                  end if
c$$$  if( (grpi.ne.grpj).OR.((j.lt.jlow).OR.(j.gt.jhigh))) then
c$$$  dx = h_X(i) - h_X(j)
c$$$  dy = h_Y(i) - h_Y(j)
c$$$  dz = h_Z(i) - h_Z(j)
c$$$  dr2 = dx*dx + dy*dy + dz*dz
c$$$  if (dr2 .le. Rcut2) then
c$$$  if (dr2 .gt. SMALL2) then
c$$$  dr = sqrt(dr2)
c$$$  else
c$$$  dr = SMALL
c$$$  end if
c$$$  Ec_HH = Ec_HH + k_e*h_chg(i)*
c$$$  +                 h_chg(j)*dexp(-debye*dr)/dr
c$$$  end if
c$$$  
c$$$  dr = sqrt(dr2)
c$$$  
c$$$  if(dr .le. vdw_cut) then
c$$$  revd = (revd_hh(i) + revd_hh(j))/2
c$$$  s1 = revd
c$$$  s2 = revd
c$$$  Ev_HH = Ev_HH + 
c$$$  +                  k_ex*((s1/dr)**12-(s2/dr)**6)
c$$$  endif   
c$$$  end if
cccccccccccccccccccccccc
ccc Modification: TONI (LH concentration)
               endif !!! LHboundbds(j)
            end do
         endif  !!! LHboundbds(i)
cccccccccccccccccccccccc
      end do        
cccccccccccccccccccc
      Ec_HC = 0.0
      Ec_HL1 = 0.0 ! adjacent linker histone/linker DNA energy
      Ec_HL2 = 0.0 ! nonadjacent linker histone/linker DNA energy
c     linker histone - DNA linker/nucleosome core interactions
      do i = h_par_start, h_par_stop
C     do i = 1,h_n  
         if(LHboundbds(i)) then
            do j = 1, n	  	  
               j1 = 3*(j-1) + 1
               j2 = j1 + 1
               j3 = j2 + 1
               dx = h_X(i) - r(j1)
               dy = h_Y(i) - r(j2)
               dz = h_Z(i) - r(j3)
               dr2 = dx*dx + dy*dy + dz*dz
               dr = sqrt(dr2)
               
               if ( dr .le. Rcut+6.0 ) then ! dr<rcut+6 i.e rcut+core radius
                  if (dr2 .gt. SMALL2) then
                     dr = sqrt(dr2)
                  else
                     dr = SMALL
                  end if
c     ... Linker DNA (type = 0)
                    if ( type(j) .eq. 0 ) then
                            if (dr. le. Rcut) then
ccccccccccccccccc
ccc   Modification: TONI (LHref)
ccc   Adapt for variable LHs the parent/non-parent condition

                                if (((i-1)/nbLH.eq.indexj(j)-1).or.
     +                          ((i-1)/nbLH-1.eq.indexj(j)-1)) then

                                Ec_HL1 = Ec_HL1 + k_e * 
     +                          h_chg(i)*q_l*dexp(-debye*dr)/dr
                                else

                                   if(fnonparLH) then       
                                      Ec_HL2 = Ec_HL2 + k_e * 
     +                                h_chg(i)*q_l*dexp(-debye*dr)/dr		       
                                   end if

                                endif
                             endif
                             if (dr .le. vdw_cut) then
ccccccccccccccc
ccc   Modification: TONI (LHref)
ccc   Adapt to variable LHs the exclude volume
                                revd = revd_hl(i)
                                s1 = revd
                                s2 = revd                  


cccccccccccccccccccc
ccc   Modification: TONI
ccccccccccccccccc
ccc   Modification: TONI (LHref)
ccc   Adapt for variable LHs the parent/non-parent condition
ccc   Prefactor for the LH/linker DNA Ev
                                    if (((i-1)/nbLH.eq.indexj(j)-1).or.
     +                               ((i-1)/nbLH-1.eq.indexj(j)-1)) then
                                    else
                                      Ev_HL2 = Ev_HL2 + 
     +                                k_exhl*((s1/dr)**12-(s2/dr)**6)
                                    end if

                               end if
                    else

                            if ((i-1)/nbLH.ne.indexj(j)-1) then

                                do k = 1,Nq
                                           k1 = 3*(k-1) + 1
                                           k2 = k1 + 1
                                           k3 = k2 + 1
                             dx = h_X(i) - ( r(j1) + a(j1)*core_pos(k1)
     +                                          + b(j1)*core_pos(k2)
     +                                          + c(j1)*core_pos(k3) )
                             dy = h_Y(i) - ( r(j2) + a(j2)*core_pos(k1)
     +                                          + b(j2)*core_pos(k2)
     +                                          + c(j2)*core_pos(k3) )
                             dz = h_Z(i) - ( r(j3) + a(j3)*core_pos(k1)
     +                                          + b(j3)*core_pos(k2)
     +                                          + c(j3)*core_pos(k3) )
                                           dr2 = dx*dx + dy*dy + dz*dz
                                           if (dr2 .gt. SMALL2) then
                                              dr = sqrt(dr2)
                                           else
                                              dr = SMALL
                                           end if

                                           if ( dr .le. Rcut ) then
!     To exclude intentionally elec. interactions between LH and nucleosome 
!     pseudo-charges, we don't add Ec_HC to the total energy:         
ccc   Modification: TONI
c     We use a flag to control if this interaction is active or inactive
                                              if(fnonparLH) then       
                                           Ec_HC = Ec_HC + k_e*h_chg(i)* 
     +                                      core_q(k)*dexp(-debye*dr)/dr
                                              end if
ccc   Old commented code
c$$$  Ec_HC = Ec_HC + k_e*h_chg(i)* 
c$$$  +                       core_q(k)*dexp(-debye*dr)/dr
                                           endif
ccccccccccccc
                                           if ( dr .le. vdw_cut ) then
ccccccccccccccc
ccc   Modification: TONI (LHref)
ccc   Adapt to variable LHs the exclude volume
                                              revd = revd_hc(i)
                                              s1 = revd
                                              s2 = revd                  

ccc   if (mod(i,nbLH).eq.1) then ! inner LH bead
ccc   Old code
c$$$  if (mod(i,3).eq.1) then ! inner LH bead
c$$$  s1 = evd_hcG
c$$$  s2 = evd_hcG	     
c$$$  else ! outer LH bead
c$$$  s1 = evd_hcC
c$$$  s2 = evd_hcC
c$$$  endif
ccccccccccccccc

cccccccccccccccccccc
ccc   Modification: TONI
                                              Ev_HC = Ev_HC + 
     +                                    k_ex*((s1/dr)**12-(s2/dr)**6)
c     Old
ccc   Ev = Ev + k_ex*((s1/dr)**12-(s2/dr)**6)
c     Test
c     e_tmp = k_ex * ((s1/dr)**12-(s2/dr)**6)
c     write(*,*)'Ev_HC',e_tmp
cccccccccccccccccccc
                                           end if
                                        end do
                                     endif
                                  end if
                               end if ! dr<Rcut+11
                            end do
cccccccccccccccccccccccc
ccc Modification: TONI (LH concentration)
                             endif  !!! LHboundbds(i)

                             
                             
                             
cccccccccccccccccccccccc
      end do

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     linker histone - tail interactions
      Ec_TH = 0.0
      Ev_HT = 0.0
      do i = h_par_start, h_par_stop, 1
C     do i = 1,h_n  
         if(LHboundbds(i)) then
            do j = 1, t_n	  	  
               dx = h_X(i) - t_X(j)
               dy = h_Y(i) - t_Y(j)
               dz = h_Z(i) - t_Z(j)
               dr2 = dx*dx + dy*dy + dz*dz
               dr = sqrt(dr2)
               if ( dr .le. Rcut) then !
                  if (dr2 .gt. SMALL2) then
                     dr = sqrt(dr2)
                  else
                     dr = SMALL
                  end if
                  if (dr. le. Rcut) then
                     Ec_TH = Ec_TH + k_e * 
     +                    h_chg(i)*t_chg(j)*dexp(-debye*dr)/dr
                  endif
               end if 
ccccccccccccccc
ccc   Modification: TONI (LHref)
ccc   Depending on the resolution, LH can have negative patches
               if ( dr .le. vdw_cut ) then
                  revd = revd_ht(i)
                  s1 = revd
                  s2 = revd                  
                  Ev_HT = Ev_HT + 
     +                 k_ex*((s1/dr)**12-(s2/dr)**6)
c     Test
c     e_tmp = k_ex * ((s1/dr)**12-(s2/dr)**6)
c     write(*,*)'Ev_HT',e_tmp
               end if
cccccccccccccccccccc
            end do
cccccccccccccccccccccccc
ccc Modification: TONI (LH concentration)
         endif  !!! LHboundbds(i)
cccccccccccccccccccccccc
      end do
      
      if ((n_c.eq.1).and.(n.eq.1)) then !<<<<<<<<
         Ev = 0.0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	 Ec_CC = 0.0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<
	 Ec_LC = 0.0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<
	 Ec_LL = 0.0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<
	 ssle = 0.0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	 ssbe = 0.0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	 sstw = 0.0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      endif !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      
C      end if   !..............................................END WITHLINK IF
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Part 3: Tail Model Internal Energy: bond + angle
      t_Eb = 0.0d0
      t_Ea = 0.0d0
c ... t_Eb = 1/2 * kb * (bi - b0)^2
      do i = t_par_start, t_par_stop - 1
C     do i = 1, t_n - 1   
        if ( t_grp(i+1) .eq. t_grp(i) ) then
          dx = t_X(i) - t_X(i+1)
          dy = t_Y(i) - t_Y(i+1)
          dz = t_Z(i) - t_Z(i+1)
          dr2 = dx*dx + dy*dy + dz*dz
          t_bond(i) = SMALL
          if (dr2 .gt. SMALL2) t_bond(i) = sqrt( dr2 )
          t_Eb = t_Eb +
     +           0.5*t_bond_c(i)*(t_bond(i) - t_bond_v(i))**2
        end if
      end do

c ... t_Ea = 1/2 * ka * (ai - a0)^2
      do i = t_par_start, t_par_stop - 2
C     do i = 1, t_n - 2  
c         to avoid calculations on non-sense inter-group angles
          if ( ( t_grp(i+1) .eq. t_grp(i) ) .and. 
     +         ( t_grp(i+2) .eq. t_grp(i) ) ) then
              x10 = t_X(i) - t_X(i+1)
              y10 = t_Y(i) - t_Y(i+1)
              z10 = t_Z(i) - t_Z(i+1)
              r10s = x10*x10 + y10*y10 + z10*z10
              if (r10s .gt. SMALL2)  then 
                  r10 = sqrt(r10s)
              else
                  r10 = SMALL
              end if

              x12 = t_X(i+2) - t_X(i+1)
              y12 = t_Y(i+2) - t_Y(i+1)
              z12 = t_Z(i+2) - t_Z(i+1)
              r12s = x12*x12 + y12*y12 + z12*z12
              if (r12s .gt. SMALL2)  then 
                  r12 = sqrt(r12s)
              else
                  r12 = SMALL
              end if

              phi = ( x10*x12 + y10*y12 + z10*z12 ) / ( r10*r12 )
              t_angle(i) = acos( phi )
              t_Ea = t_Ea + 0.5d0 * t_angle_c(i) * 
     +               ( t_angle(i) - t_angle_v(i) )**2
          end if
      end do 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Part 4: Tail Model Internal Energy: electrostatics + excluded volume
c     A.  On different tails
c     B.  On the same tail: j > i+2 (consistent with UHBD derivation)
c     C.  But fixed tails on the same core should not have interactions
c             (not a problem after using stretching between fixed tails 
c              and their own core)
      do i = t_par_start, t_par_stop - 1
C     do i = 1, t_n - 1  
          do j = i + 1, t_n
              if (  ( t_grp(i) .ne. t_grp(j) ) 
     +              .or. ( j .gt. (i+2) )      ) then

                  dx = t_X(i) - t_X(j)
                  dy = t_Y(i) - t_Y(j)
                  dz = t_Z(i) - t_Z(j)
                  dr2 = dx*dx + dy*dy + dz*dz

c ... t_Ec = k_e * qi * qj * dexp(-debye*dr) / dr
                  if (dr2 .le. Rcut2) then
                      if (dr2 .gt. SMALL2) then
                          dr = sqrt(dr2)
                      else
                          dr = SMALL
                      end if
                      if ( (i-1)/ntpc .eq. (j-1)/ntpc ) then
                          Ec_TT1 = Ec_TT1 + k_e*t_chg(i)*t_chg(j)*
     +                             dexp(-debye*dr) / dr
cccccccccccccccccc
ccc Modification: DURBA and TONI
c     This interaction was missing in the force-field of the total energy,
c     but was present in the 'partialpot.f', i.e., it was evaluated during
c     the Monte Carlo cycles. Acordingly, we added it also here.
                      else
                         Ec_TT2 = Ec_TT2 + k_e*t_chg(i)*t_chg(j)*
     +                             dexp(-debye*dr) / dr
cccccccccccccccccc
                      end if
                  end if

c ... t_Ev = e(i,j) * [3(r0(i,j)/r(i,j))^8 - 4(r0(i,j)/r(i,j))^6 + 1]
                  if (dr2 .le. vdw_cut2) then
                      if (dr2 .gt. SMALL2) then
                          dr = sqrt(dr2)
                      else
                          dr = SMALL
                      end if
cccccccccccccccccccc
ccc Modification: TONI
                      if ( (i-1)/ntpc .eq. (j-1)/ntpc ) then
                        Ev_TT1 = Ev_TT1 + 
     +                       k_ex*((t_exv_d/dr)**12-(t_exv_d/dr)**6)
                      else
                        Ev_TT2 = Ev_TT2 + 
     +                       k_ex*((t_exv_d/dr)**12-(t_exv_d/dr)**6)
                      end if
c     Old
c                      t_Ev = t_Ev + k_ex * 
c     +                   ((t_exv_d/dr)**12-(t_exv_d/dr)**6)
c     Test
c                      e_tmp = k_ex*((t_exv_d/dr)**12-(t_exv_d/dr)**6)
c                      write(*,*)'Ev_TT',e_tmp
cccccccccccccccccccc
                  end if

 1000         end if
          end do
      end do

c ... adding to total Ec and Ev
      Ec_TT = Ec_TT1 + Ec_TT2
      t_Ec  = Ec_TT
ccccccccccccccccccc
ccc Modification: TONI
      Ev_TT = Ev_TT1 + Ev_TT2
      t_Ev = Ev_TT
c     Ev = Ev + t_Ev
ccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Part 5: Tail External Interaction Energy with Core and Linker DNA: 
c         electrostatics + excluded volume
c             + stretching (fixed tails with their own core)
      do i = t_par_start, t_par_stop
          
          do j = 1, n
              j1 = 3*(j-1) + 1
              j2 = j1 + 1
              j3 = j2 + 1
              dx = t_X(i) - r(j1)
              dy = t_Y(i) - r(j2)
              dz = t_Z(i) - r(j3)
              dr2 = dx*dx + dy*dy + dz*dz
	      dr = sqrt(dr2)

              if ( dr .le. Rcut+6.0 ) then ! dr<rcut+6 i.e rcut+core radius
                  if (dr2 .gt. SMALL2) then
                      dr = sqrt(dr2)
                  else
                      dr = SMALL
                  end if
c ... Linker DNA (type = 0)
                  if ( type(j) .eq. 0 ) then
		      if (dr. le. Rcut) then
                          Ec_TL = Ec_TL + k_e * 
     +                       t_chg(i)*q_l*dexp(-debye*dr)/dr
     		      endif
                      if (dr .le. vdw_cut) then
                          s1 = evd_tl
                          s2 = evd_tl
cccccccccccccccccccc
ccc Modification: TONI
                          Ev_TL = Ev_TL + 
     +                         k_ex*((s1/dr)**12-(s2/dr)**6)
c     Old
ccc                   Ev = Ev + k_ex*((s1/dr)**12-(s2/dr)**6)
c     Test
c                          e_tmp = k_ex * ((s1/dr)**12-(s2/dr)**6)
c                          write(*,*)'Ev_TL',e_tmp
cccccccccccccccccccc
                      end if
                  else
c=============Core: start================c  
c ... Core (type = 1 or 2; 2 for last core)
                     if (( t_fix(i) .eq. 1 ) .and. 
     + 			((i-1)/ntpc .eq. indexj(j)-1 )) then
c         Fixed tail beads only have stretching with their own core.
c         find x0, y0, z0 - the coordinates where the fixed tail should be.
          		 x0 = r(j1) + a(j1) * t_X0(i)
     +               		+ b(j1) * t_Y0(i)
     +               		+ c(j1) * t_Z0(i)
          		 y0 = r(j2) + a(j2) * t_X0(i)
     +               		+ b(j2) * t_Y0(i)
     +               		+ c(j2) * t_Z0(i)
          		 z0 = r(j3) + a(j3) * t_X0(i)
     +               		+ b(j3) * t_Y0(i)
     +               		+ c(j3) * t_Z0(i)
c         compute stretching energy
          		dx = x0 - t_X(i)
          		dy = y0 - t_Y(i)
          		dz = z0 - t_Z(i)
ccccccccccccccccccccc
ccc Modification: DURBA and TONI
c     The prefactor 0.5 has to be double precision (compatibility with master01)
c          		Eb_tc = Eb_tc + 0.5*h_tc*(dx*dx+dy*dy+dz*dz)
          		Eb_tc = Eb_tc + 0.50d0*h_tc*(dx*dx+dy*dy+dz*dz)
ccccccccccccccccccccc
          		go to 8000
      		    
                  else
          		do k = 1,Nq
              		  k1 = 3*(k-1) + 1
              		  k2 = k1 + 1
              		  k3 = k2 + 1
              		  dx = t_X(i) - ( r(j1) + a(j1)*core_pos(k1)
     +                              + b(j1)*core_pos(k2)
     +                              + c(j1)*core_pos(k3) )
              		  dy = t_Y(i) - ( r(j2) + a(j2)*core_pos(k1)
     +                              + b(j2)*core_pos(k2)
     +                              + c(j2)*core_pos(k3) )
              		  dz = t_Z(i) - ( r(j3) + a(j3)*core_pos(k1)
     +                              + b(j3)*core_pos(k2)
     +                              + c(j3)*core_pos(k3) )
              		  dr2 = dx*dx + dy*dy + dz*dz
              		  if (dr2 .gt. SMALL2) then
                  		dr = sqrt(dr2)
              		  else
                  		dr = SMALL
              		  end if
			    if (dr.le.Rcut) then
                            if ( (i-1)/ntpc .eq. indexj(j)-1 ) then
                               Ec_CT1 = Ec_CT1 + k_e * 
     +                         t_chg(i)*core_q(k)*dexp(-debye*dr)/dr
                            else !<<<<<<<<<<<<comment below for Ect2=0
                               Ec_CT2 = Ec_CT2 + k_e * 
     +                         t_chg(i)*core_q(k)*dexp(-debye*dr)/dr
                            end if
			    endif
                          if ( dr .le. vdw_cut ) then
                            s1 = evd_tc
                            s2 = evd_tc
cccccccccccccccccccc
ccc Modification: TONI
                            if ( (i-1)/ntpc .eq. (j-1)/nt ) then
                               Ev_CT1 = Ev_CT1 + 
     +                              k_ex*((s1/dr)**12-(s2/dr)**6)
                             else !<<<<<<<<<<<<comment below for Ect2=0
                               Ev_CT2 = Ev_CT2 + 
     +                              k_ex*((s1/dr)**12-(s2/dr)**6)
                             end if
c     Old
ccc                   Ev = Ev + k_ex*((s1/dr)**12-(s2/dr)**6)
c     Test
c                            e_tmp = k_ex * ((s1/dr)**12-(s2/dr)**6)
c                            write(*,*)'Ev_CT',e_tmp,e_tmp1,e_tmp2
cccccccccccccccccccc
                          end if
                        end do
                    end if
c=============Core: end===================c
 8000             end if
              end if  ! dr<Rcut+5
          end do
      end do
 
      Ec_CT = Ec_CT1 + Ec_CT2
      Ev_CT = Ev_CT1 + Ev_CT2
      Ev = Ev_LL + Ev_CC + Ev_TT + Ev_LC + Ev_CT 
     +     + Ev_TL + Ev_HL1 + Ev_HL2 + Ev_HC + Ev_HT
     +     + Ev_HH

c     The flag controls if we add the non-parental interactions
      if(fnonparLH) then
         Ec = Ec_LL + Ec_CC + Ec_TT + Ec_LC + Ec_CT 
     +        + Ec_TL + Ec_HL1 + Ec_HL2 + Ec_HC + Ec_TH
     +        + Ec_HH
      else
         Ec = Ec_LL + Ec_CC + Ec_TT + Ec_LC + Ec_CT + Ec_TL   
     +        + Ec_HL1 +  Ec_TH + Ec_HH
      end if
      
        pot_send_array(1)=Ev_LL
        pot_send_array(2)=Ev_CC
        pot_send_array(3)=Ev_TT
        pot_send_array(4)=Ev_LC
        pot_send_array(5)=Ev_CT
        pot_send_array(6)=Ev_TL
        pot_send_array(7)=Ev_HL1
        pot_send_array(8)=Ev_HL2
        pot_send_array(9)=Ev_HC
        pot_send_array(10)=Ev_HT
        pot_send_array(11)=Ev_HH
        pot_send_array(12)=Ec_LL
        pot_send_array(13)=Ec_CC
        pot_send_array(14)=Ec_TT
        pot_send_array(15)=Ec_LC
        pot_send_array(16)=Ec_CT
        pot_send_array(17)=Ec_TL
        pot_send_array(18)=Ec_HL1
        pot_send_array(19)=Ec_HL2
        pot_send_array(20)=Ec_HC
        pot_send_array(21)=Ec_TH
        pot_send_array(22)=Ec_HH
        pot_send_array(23)=t_Eb
        pot_send_array(24)=t_Ea
        pot_send_array(25)=Eb_tc
        pot_send_array(26)=Ec
        pot_send_array(27)=Ev
        pot_send_array(28)=ELHel

      call MPI_ALLREDUCE(pot_send_array,sum_pot_send_array,28,
     +MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        Ev_LL=sum_pot_send_array(1)
        Ev_CC=sum_pot_send_array(2)
        Ev_TT=sum_pot_send_array(3)
        Ev_LC=sum_pot_send_array(4)
        Ev_CT=sum_pot_send_array(5)
        Ev_TL=sum_pot_send_array(6)
        Ev_HL1=sum_pot_send_array(7)
        Ev_HL2=sum_pot_send_array(8)
        Ev_HC=sum_pot_send_array(9)
        Ev_HT=sum_pot_send_array(10)
        Ev_HH=sum_pot_send_array(11)
        Ec_LL=sum_pot_send_array(12)
        Ec_CC=sum_pot_send_array(13)
        Ec_TT=sum_pot_send_array(14)
        Ec_LC=sum_pot_send_array(15)
        Ec_CT=sum_pot_send_array(16)
        Ec_TL=sum_pot_send_array(17)
        Ec_HL1=sum_pot_send_array(18)
        Ec_HL2=sum_pot_send_array(19)
        Ec_HC=sum_pot_send_array(20)
        Ec_TH=sum_pot_send_array(21)
        Ec_HH=sum_pot_send_array(22)
        t_Eb=sum_pot_send_array(23)
        t_Ea=sum_pot_send_array(24)
        Eb_tc=sum_pot_send_array(25)
        Ec=sum_pot_send_array(26)
        Ev=sum_pot_send_array(27)
        ELHel=sum_pot_send_array(28)
                      
       

C Sum up all the energy terms
      E(1) = ssle
      E(2) = ssbe
      E(3) = sstw
      E(4) = Ec
      E(5) = Ev
      E(6) = E(1)+E(2)+E(3)+E(4)+E(5)+t_Eb+t_Ea+Eb_tc+ELHel

       close(unit=133)
 1089  format(1x,A45)

 1809 format(1x, A25, 2x,F10.6)

 2089 format(1x,A15,I2,2x,A15,I2) 
 2289 format(1x,A12,I3,2x,A4,2x,F10.6)
      
      
      
         end


