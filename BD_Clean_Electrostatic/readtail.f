      subroutine readtail(tailfile,n_c,n, t_n, 
     +             t_grp, t_fix, t_chg, t_rad, t_mass,
     +             t_X, t_Y, t_Z,
     +             t_X0, t_Y0, t_Z0,
     +             t_bond, t_bond_v, t_bond_c, 
     +             t_angle, t_angle_v, t_angle_c,n3, r, a, b, c,  
     +             h_n, h_X, h_Y, h_Z, h_chg,ro, qH1G, qH1C1, qH1C2,
     +             h_X0,h_Y0,h_Z0,withlink,myid,nlb,type)

ccccccccccccccccccccc
ccc Modification: TONI (LHref)
ccc   We make use of global parameters/variables
      ! nbNt,nbGh,nbCt,nbLH
      ! rxyzNt,rxyzGh,rxyzCt: LH coordinates
      ! rqNt,rqGh,rqCt: LH charges
      ! revd_hh*,revd_hc*,revd_hl*: LH excluded volume
      ! rkstr*,rkben*,rktor*: LH elastic parameters
      use modglob
ccccccccccccccccccccc
 
      implicit NONE

c n_c: number of cores, t_n: total tail segments, nt: number of tail segments per core
c N_unit: number of beads of core and linker DNA per unit (monomer)
      integer n_c, t_n, h_n,n, N_unit,myid,type(n)
      character*10 resname(t_n)
      integer t_grp(t_n), t_fix(t_n), nlb(n_c)
      double precision t_chg(t_n), t_rad(t_n), t_mass(t_n), ro
      double precision t_X(t_n), t_Y(t_n), t_Z(t_n)
      double precision t_X0(t_n), t_Y0(t_n), t_Z0(t_n)
      double precision t_bond(t_n), t_bond_v(t_n), t_bond_c(t_n)
      double precision t_angle(t_n), t_angle_v(t_n), t_angle_c(t_n)
      double precision h_X(h_n),h_Y(h_n),h_Z(h_n),h_chg(h_n)
      double precision h_X0(3),h_Y0(3),h_Z0(3)
      integer n3, count
      double precision r(n3), a(n3), b(n3), c(n3),dist
      character   tailfile*50

c ... ntgrp is the amout of tails (groups) on each core
c ... gstart, gend, and gfix are the starting bead, the ending bead, 
c ...     and the fixed bead of each tail, respectively. 
      integer i, j, ntpc, k, ki, m1, m2, m3, ntgrp, gstart, gend, gfix
      double precision dx, dy, dz, dr2, dr, dbx, dby, dbz, dbr
      double precision x12, y12, z12, x32, y32, z32, p123, r12, r32
      double precision pi, xrel, yrel, zrel, qH1G, qH1C1, qH1C2
      parameter(pi = 3.14159265358979d0)

ccccccccc
ccc Modification: TONI (LHref)
      double precision rq
      integer index,index2
      double precision rehh,rehc,rehl,reht
      double precision kstr,kben,ktor
      integer j1,j2
      double precision x1,x2,y1,y2,z1,z2
      double precision rxyz1(3),rxyz2(3),rxyz3(3),drxyz(3)
      double precision rxyz12(3),rxyz23(3)
      double precision norm12,norm23,dotpr
      double precision cosbet,beta_tmp
      integer k0,k1,k2,k3
ccccccccc
      logical withlink
     
      ntpc = t_n/n_c
      N_unit = n3/3/n_c
      
!      qH1G = 12.4 ! 200mM ! 0.0!  4.6 ! 010mM  
!      qH1C = 59.8 ! 200mM ! 0.0 ! 16.7 ! 010mM 

c Load tail model data file      
      open (unit = 1, file = tailfile, form = 'formatted',
     +      access = 'sequential', status = 'old')
      write(*,'(A,i2)'),' TAIL FILE BEING READ BY.......',myid

      do 50 i = 1,ntpc
        read (unit = 1, fmt = 40, end = 100)
     +        resname(i), 
     +        t_grp(i), t_fix(i), t_chg(i), t_rad(i), t_mass(i),
     +        t_X0(i), t_Y0(i), t_Z0(i),
     +        t_bond_v(i), t_bond_c(i),
     +        t_angle_v(i), t_angle_c(i)
   40   format( A10, 
     +          I10, I10, F10.2, F10.2, F10.2,
     +          F10.3, F10.3, F10.3,
     +          F10.2, F10.2,
     +          F10.2, F10.2 )
   50 continue

      close ( unit = 1 )


c Unit conversions
      do i = 1,ntpc
c ....... Radius: A => nm
          t_rad(i) = t_rad(i) / 10.0d0
c ....... Mass: g/mol => (*** mass is not used in BD ***)
c
c ....... X, Y, Z: A => nm
          t_X0(i) = t_X0(i) / 10.0d0
          t_Y0(i) = t_Y0(i) / 10.0d0
          t_Z0(i) = t_Z0(i) / 10.0d0
c ....... Bond value: A => nm
          t_bond_v(i) = t_bond_v(i) / 10.0d0
c ....... Bond constant: kcal/mol/A^2 => kcal/mol/nm^2
          t_bond_c(i) = t_bond_c(i) * 100.0d0
c ....... Angle value: Deg => rad
          t_angle_v(i) = t_angle_v(i) * pi / 180.0d0
      end do

c Set initial bond lengths
      do i = 1, ntpc - 1
        if ( t_grp(i+1) .eq. t_grp(i) ) then
          dx = t_X0(i) - t_X0(i+1)
          dy = t_Y0(i) - t_Y0(i+1)
          dz = t_Z0(i) - t_Z0(i+1)
          dr2 = dx*dx + dy*dy + dz*dz
          t_bond(i) = 0.0d0
          if (dr2 .ne. 0.0d0) t_bond(i) = sqrt( dr2 )
        else
          t_bond(i) = 0.0d0
        end if
      end do
      t_bond(ntpc) = 0.0d0
      
c Set initial angle values
c ... acos() return unit of rad, so no unit conversion is needed
      do i = 1,ntpc-2
         x12 = t_X0(i) - t_X0(i+1)
         y12 = t_Y0(i) - t_Y0(i+1)
         z12 = t_Z0(i) - t_Z0(i+1)
         r12 = x12*x12 + y12*y12 + z12*z12
         if (r12 .ne. 0.0d0)  r12 = sqrt(r12)

         x32 = t_X0(i+2) - t_X0(i+1)
         y32 = t_Y0(i+2) - t_Y0(i+1)
         z32 = t_Z0(i+2) - t_Z0(i+1)
         r32 = x32*x32 + y32*y32 + z32*z32
         if (r32 .ne. 0.0d0)  r32 = sqrt(r32)

         p123 = ( x12*x32 + y12*y32 + z12*z32 ) / ( r12*r32 )
         t_angle(i) = acos( p123 )
      enddo 
      t_angle(ntpc-1) = 0.0d0
      t_angle(ntpc) = 0.0d0


c Set rest of parameters in ntpc to t_n
c ntgrp is the number of tails (tail groups) on each core
      ntgrp = t_grp(ntpc)
      do j = 2, n_c
          k = (j-1) * ntpc
          do i = 1, ntpc
              ki = k + i
              resname(ki) = resname(i)
              t_grp(ki) = t_grp(i) + (j-1)*ntgrp
              t_fix(ki) = t_fix(i)
              t_chg(ki) = t_chg(i)
              t_rad(ki) = t_rad(i)
              t_mass(ki) = t_mass(i)

              t_X0(ki) = t_X0(i)
              t_Y0(ki) = t_Y0(i)
              t_Z0(ki) = t_Z0(i)

              t_bond(ki) = t_bond(i)
              t_bond_v(ki) = t_bond_v(i)
              t_bond_c(ki) = t_bond_c(i)
              t_angle(ki) = t_angle(i)
              t_angle_v(ki) = t_angle_v(i)
              t_angle_c(ki) = t_angle_c(i)
          end do
      end do

c Transform t_X, t_Y, and t_Z with r, a, b, and c
      index=1
      do j = 1, n_c
          k = (j-1) * ntpc
          m1 = (index-1) * 3 + 1
          m2 = m1 + 1
          m3 = m2 + 1
          do i = 1, ntpc
              ki = k + i
              t_X(ki) = r(m1) + a(m1) * t_X0(ki)
     +                        + b(m1) * t_Y0(ki)
     +                        + c(m1) * t_Z0(ki)
              t_Y(ki) = r(m2) + a(m2) * t_X0(ki)
     +                        + b(m2) * t_Y0(ki)
     +                        + c(m2) * t_Z0(ki)
              t_Z(ki) = r(m3) + a(m3) * t_X0(ki)
     +                        + b(m3) * t_Y0(ki)
     +                        + c(m3) * t_Z0(ki)
          end do
          index=index+nlb(j)+1         
      end do
      



      if (withlink) then
ccc Setup the LH positions, charges, and VderW
ccc Notice: now the initial relative positions are stored in rxyzNt,rxyzGh,rxyzCt
ccc !!! We initially generate the same h_X vector... but we'll need vectors
ccc for the different domains to introduce flexibility
      index = 0
      index2 = 1
      do j = 1, n_c
         m1 = (index2-1)*3+1
         m2 = m1+1
         m3 = m2+1
ccc LH
         do k = 1,nbLH
            index = index + 1
            xrel = rxyzLH(k,1)
            yrel = rxyzLH(k,2)
            zrel = rxyzLH(k,3)
            h_X(index) = r(m1) + a(m1)*xrel + b(m1)*yrel + c(m1)*zrel
            h_Y(index) = r(m2) + a(m2)*xrel + b(m2)*yrel + c(m2)*zrel
            h_Z(index) = r(m3) + a(m3)*xrel + b(m3)*yrel + c(m3)*zrel
            
            rq = rqLH(k)
            h_chg(index) = rq

            rehh = revd_hhLH(k)
            rehc = revd_hcLH(k)
            rehl = revd_hlLH(k)
            reht = revd_htLH(k)
            revd_hh(index) = rehh
            revd_hc(index) = rehc
            revd_hl(index) = rehl
            revd_ht(index) = reht
            
            kstr = rkstrLH(k)
            kben = rkbenLH(k)
            ktor = rktorLH(k)
            rkstr(index) = kstr
            rkben(index) = kben
            rktor(index) = ktor

            h_conn(index) = connLH(k)
            h_bend(index) = bendLH(k)
         enddo
         index2=index2+nlb(j)+1         
      enddo

ccc For numerical reasons we store the elastic constants divided by two
      rkstrd2 = rkstr/2.
      rkbend2 = rkben/2.
      rktord2 = rktor/2.



ccc Computationally it's convenient to store the equilibrium positions 
ccc in common vectors for all LH beads
      do i = 1,h_n
         j = modulo(i,nbLH)
         if (j.EQ.0) j = nbLH
         !!! Streching 
         rkstreq(i) = rkstreq(j)
         !write(*,*) 'Strech eq',i,j,rkstreq(i)
         !!! Bending
         rkbeneq(i) = rkbeneq(j)
         !write(*,*) 'Ben_eq',i,j,rkbeneq(i)        
      enddo

ccc Check
c$$$      do j =1,nbLH
c$$$         xrel = h_X(j)
c$$$         yrel = h_Y(j)
c$$$         zrel = h_Z(j)
c$$$         rq = h_chg(j)
c$$$         rehh = revd_hh(j)
c$$$         rehc = revd_hc(j)
c$$$         rehl = revd_hl(j)
c$$$         reht = revd_ht(j)
c$$$         write(*,*) xrel,yrel,zrel,rq,rehh,rehc,rehl,reht
c$$$         kstr = rkstr(j)
c$$$         kben = rkben(j)
c$$$         ktor = rktor(j)
c$$$         write(*,*) 'Elast:',kstr,kben,ktor
c$$$         kstr = rkstrd2(j)
c$$$         kben = rkbend2(j)
c$$$         ktor = rktord2(j)
c$$$         write(*,*) 'Elast:',kstr,kben,ktor
c$$$      enddo


ccc LH groups
      do i = 1,h_n
         j = (i-1)/nbLH + 1
         LH_grp(i) = j
c         write(*,*) 'LH_grp',i,LH_grp(i)
      enddo

ccc To avoid the modification of the calling functions we keep the h_X0 vectors
ccc but initialized to zero
ccc Notice that they're not used in the current code
      xrel =  0
      yrel =  0
      zrel =  0
      h_X0(1)=xrel
      h_Y0(1)=yrel
      h_Z0(1)=zrel

      xrel =  rxyzCt(1,1)
      yrel =  rxyzCt(1,2)
      zrel =  rxyzCt(1,3)
      h_X0(2)= 0
      h_Y0(2)= 0
      h_Z0(2)= 0
      
      xrel = 0
      yrel = 0
      zrel = 0
      h_X0(3)=xrel
      h_Y0(3)=yrel
      h_Z0(3)=zrel

      endif
     
c >>>>>>>>>>>>>>>>>>> set H4 tails to zero >>>>>>>>>>>>>>>>>>>>>>>>
!      do j = 1, n_c
!         k = (j-1)*ntpc
!         do i = 35, 44 !H3: 1-16, H4: 17-26, H2AN: 27-34, H2B:35-44 
!      	    ki = k+i
!	    t_chg(ki) = 0.0
!      	 enddo
!      enddo
!      write(*,*) 'WARNING WARNING !! A tail is being neutralized'
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      

      return

100   print*, 'Premature end of file in subroutine: readtail'

      end







