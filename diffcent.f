Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
      subroutine   diffcent(n_c,nc3,n,n3,type,r,ro,d1,go,a,b,c,
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
     +force_dnaStr_X_num,force_dnaStr_Y_num,force_dnaStr_Z_num,
     +force_dnaBend_X_num,force_dnaBend_Y_num,force_dnaBend_Z_num,
     +force_dnaTwist_X_num,force_dnaTwist_Y_num,force_dnaTwist_Z_num,
     +                force_dnaStr_X,force_dnaStr_Y,force_dnaStr_Z,  !analytic force values inside the function
     +                force_dnaBend_X,force_dnaBend_Y,force_dnaBend_Z,
     +             force_dnaTwist_X,force_dnaTwist_Y,force_dnaTwist_Z,
     +np,myid,ierr)
      
      ! Numerical Partial differentiation of energy function
      ! perturbation is done only on a specified coordiate but for all beads

      implicit none  !CSUN 6/2/04
c        Variables needed for the update_mod routine
      integer n_c,nc3,n,num,n3, type(n),myid,ierr,np,index,j1,j2,j3,h_n
      double precision r(n3),ro,d1,go, a(n3),b(n3),c(n3)
      double precision aold(n3),bold(n3),cold(n3)
      double precision alpha(n),beta(n),gamma(n), length(n)
      double precision a_dna(nc3),b_dna(nc3),c_dna(nc3)
      double precision alpha_p(n_c),beta_p(n_c),gamma_p(n_c)

c        Variables needed for the potential calculation routine
      
      double precision lo,hd2,gd2,sd2,k_e,k_ex
      integer Nq, Nq3
      double precision core_pos(Nq3), core_q(Nq)
      double precision debye, q_l, phi_o(n_c)
      double precision E(6)
      
      integer nm1
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
      
      integer indexj(n)

        
c   Variables needed for the function output        
        
      double precision    force_dnaStr_X_num(n),force_dnaStr_Y_num(n)
      double precision force_dnaStr_Z_num(n)
      double precision force_dnaBend_X_num(n),force_dnaBend_Y_num(n)
      double precision force_dnaBend_Z_num(n)
      double precision force_dnaTwist_X_num(n),force_dnaTwist_Y_num(n)
      double precision force_dnaTwist_Z_num(n)
      
c     local variables needed by the function 
      integer i,i1,i2,i3
      double precision rs(n3),as(n3),bs(n3),cs(n3)
      !double precision rsf(n3),rsb(n3)
      double precision alphas(n),betas(n),gammas(n),lengths(n)
      double precision a_dnas(nc3),b_dnas(nc3),c_dnas(nc3)
      double precision alpha_ps(n_c),beta_ps(n_c),gamma_ps(n_c)
      double precision Es(6),Esf(6),Esb(6) !Shifted Energy backword and forward
c     local variables needed for testgh part
      double precision fc, yhy, xc(n3),gc(n3),y(n3),vec(n3)
      integer isAnalytic
!analytic DNA forces passed to compare
      double precision force_dnaStr_X(n),force_dnaStr_Y(n)
      double precision force_dnaStr_Z(n)
      double precision force_dnaBend_X(n),force_dnaBend_Y(n)
      double precision force_dnaBend_Z(n)
      double precision force_dnaTwist_X(n),force_dnaTwist_Y(n)
      double precision force_dnaTwist_Z(n)
 
      
      do 130 i = 1,n
          
        i1 = 3*(i-1) + 1
        i2 = i1 + 1
        i3 = i2 + 1
        
      ! Differentation with respect to X axis  
      ! forward step
        rs=r;
        rs(i1)=r(i1)+0.0050d0
      as=aold
      bs=bold
      cs=cold
        ! calculate the new position vectors,euler angles
      call update_mod(n_c,nc3,n,n3,type,rs,ro,d1,go,as,bs,cs,
     +        alphas,betas,gammas,lengths,a_dnas,b_dnas,c_dnas,
     +        alpha_ps,beta_ps,gamma_ps,i1+1)
     
      
      call potential(n_c,n,n3,type,rs,as,bs,cs,alphas,betas,gammas,
     +	  lengths,beta_ps,lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,Esf,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    np,myid,ierr,i1+1)
      ! backward step
        rs=r;
        rs(i1)=r(i1)-0.0050d0
      as=aold
      bs=bold
      cs=cold
        ! calculate the new position vectors,euler angles
      call update_mod(n_c,nc3,n,n3,type,rs,ro,d1,go,as,bs,cs,
     +        alphas,betas,gammas,lengths,a_dnas,b_dnas,c_dnas,
     +        alpha_ps,beta_ps,gamma_ps,i1+1)
     
      
      call potential(n_c,n,n3,type,rs,as,bs,cs,alphas,betas,gammas,
     +	  lengths,beta_ps,lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,Esb,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    np,myid,ierr,i1+1)
            
      
      force_dnaStr_X_num(i)=-(Esf(1)-Esb(1))/(2*0.0050d0)
      force_dnaBend_X_num(i)=-(Esf(2)-Esb(2))/(2*0.0050d0)
      force_dnaTwist_X_num(i)=-(Esf(3)-Esb(3))/(2*0.0050d0)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Differentation with respect to Y axis
      !forward step
        rs=r;
        rs(i2)=r(i2)+0.0050d0
      as=aold
      bs=bold
      cs=cold  
        ! calculate the new position vectors,euler angles
      call update_mod(n_c,nc3,n,n3,type,rs,ro,d1,go,as,bs,cs,
     +        alphas,betas,gammas,lengths,a_dnas,b_dnas,c_dnas,
     +        alpha_ps,beta_ps,gamma_ps,i2+1)
     
      
      call potential(n_c,n,n3,type,rs,as,bs,cs,alphas,betas,gammas,
     +	  lengths,beta_ps,lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,Esf,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    np,myid,ierr,i2+1)
      
      !backward step
        rs=r;
        rs(i2)=r(i2)-0.0050d0
      as=aold
      bs=bold
      cs=cold  
        ! calculate the new position vectors,euler angles
      call update_mod(n_c,nc3,n,n3,type,rs,ro,d1,go,as,bs,cs,
     +        alphas,betas,gammas,lengths,a_dnas,b_dnas,c_dnas,
     +        alpha_ps,beta_ps,gamma_ps,i2+1)
     
      
      call potential(n_c,n,n3,type,rs,as,bs,cs,alphas,betas,gammas,
     +	  lengths,beta_ps,lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,Esb,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    np,myid,ierr,i2+1)
      
      force_dnaStr_Y_num(i)=-(Esf(1)-Esb(1))/(2*0.0050d0)
      force_dnaBend_Y_num(i)=-(Esf(2)-Esb(2))/(2*0.0050d0)
      force_dnaTwist_Y_num(i)=-(Esf(3)-Esb(3))/(2*0.0050d0)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Differentation with respect to Z axis
      !forward step
        rs=r;
        rs(i3)=r(i3)+0.0050d0
        as=aold
      bs=bold
      cs=cold
        ! calculate the new position vectors,euler angles
      call update_mod(n_c,nc3,n,n3,type,rs,ro,d1,go,as,bs,cs,
     +        alphas,betas,gammas,lengths,a_dnas,b_dnas,c_dnas,
     +        alpha_ps,beta_ps,gamma_ps,i3+1)
     
      
      call potential(n_c,n,n3,type,rs,as,bs,cs,alphas,betas,gammas,
     +	  lengths,beta_ps,lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,Esf,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    np,myid,ierr,i3+1)
      
      !backward step
        rs=r;
        rs(i3)=r(i3)-0.0050d0
        as=aold
      bs=bold
      cs=cold
        ! calculate the new position vectors,euler angles
      call update_mod(n_c,nc3,n,n3,type,rs,ro,d1,go,as,bs,cs,
     +        alphas,betas,gammas,lengths,a_dnas,b_dnas,c_dnas,
     +        alpha_ps,beta_ps,gamma_ps,i3+1)
     
      
      call potential(n_c,n,n3,type,rs,as,bs,cs,alphas,betas,gammas,
     +	  lengths,beta_ps,lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,Esb,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    np,myid,ierr,i3+1)
      
      force_dnaStr_Z_num(i)=-(Esf(1)-Esb(1))/(2*0.0050d0)      
      force_dnaBend_Z_num(i)=-(Esf(2)-Esb(2))/(2*0.0050d0)
      force_dnaTwist_Z_num(i)=-(Esf(3)-Esb(3))/(2*0.0050d0)
          
130   continue
      
      
          open(unit=299,name='Numerical Forces.txt',access='SEQUENTIAL',
     +        status='unknown')
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!! Display forces on file
      write(unit=299, fmt=1089) 'DNA Forces:'
      write(unit=299, fmt=1089) '1-Strecthing Forces:'
      do 133 i = 1,n
          write(unit=299, fmt=1090) 'for bead #',i
          write(unit=299, fmt=1091) 'f_x = ',force_dnaStr_X_num(i),
     + 'f_y = ',force_dnaStr_Y_num(i),'f_z = ',force_dnaStr_Z_num(i)
      
133   continue
      write(unit=299, fmt=1089) '1-Bending Forces:'
      do 131 i = 1,n
          write(unit=299, fmt=1090) 'for bead #',i
          write(unit=299, fmt=1091) 'f_x = ',force_dnaBend_X_num(i),
     + 'f_y = ',force_dnaBend_Y_num(i),'f_z = ',force_dnaBend_Z_num(i)
      
131   continue
      write(unit=299, fmt=1089) '1-Twisting Forces:'
      do 132 i = 1,n
          write(unit=299, fmt=1090) 'for bead #',i
          write(unit=299, fmt=1091) 'f_x = ',force_dnaTwist_X_num(i),
     + 'f_y = ',force_dnaTwist_Y_num(i),'f_z = ',force_dnaTwist_Z_num(i)
      
132   continue
      
 
!!!!! testing the derivative using testgh routine
! test the stretching energy derivative
      xc=r
      yhy=0
      num=n
      isAnalytic=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! testing strecthing energy 
            fc=E(1)
       do 134 i = 1,n

          i1 = 3*(i-1) + 1
          i2 = i1 + 1
          i3 = i2 + 1
          y(i1)=0
          y(i2)=0
          y(i3)=0
          if (isAnalytic.eq.1) then
              gc(i1)=-force_dnaStr_X(i)
              gc(i2)=-force_dnaStr_Y(i)
              gc(i3)=-force_dnaStr_Z(i)
              else
              gc(i1)=-force_dnaStr_X_num(i)
              gc(i2)=-force_dnaStr_Y_num(i)
              gc(i3)=-force_dnaStr_Z_num(i)
           endif         
          
134    continue
       y(1)=0.05
       
      call       testgh(n3,xc,fc,gc,y,yhy,vec,
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
            fc=E(2)
       do 135 i = 1,n

          i1 = 3*(i-1) + 1
          i2 = i1 + 1
          i3 = i2 + 1
          y(i1)=0
          y(i2)=0
          y(i3)=0
          if (isAnalytic.eq.1) then
              gc(i1)=-force_dnaBend_X(i)
              gc(i2)=-force_dnaBend_Y(i)
              gc(i3)=-force_dnaBend_Z(i)
              else
              gc(i1)=-force_dnaBend_X_num(i)
              gc(i2)=-force_dnaBend_Y_num(i)
              gc(i3)=-force_dnaBend_Z_num(i)
           endif
      
                  
          
135    continue
       y(1)=0.05
      call       testgh(n3,xc,fc,gc,y,yhy,vec,
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
            fc=E(3)
       do 136 i = 1,n

          i1 = 3*(i-1) + 1
          i2 = i1 + 1
          i3 = i2 + 1
          y(i1)=0
          y(i2)=0
          y(i3)=0
          if (isAnalytic.eq.1) then
              gc(i1)=-force_dnaTwist_X(i)
              gc(i2)=-force_dnaTwist_Y(i)
              gc(i3)=-force_dnaTwist_Z(i)
              else
              gc(i1)=-force_dnaTwist_X_num(i)
              gc(i2)=-force_dnaTwist_Y_num(i)
              gc(i3)=-force_dnaTwist_Z_num(i)
           endif         
          
136    continue
       y(1)=0.05
      call       testgh(n3,xc,fc,gc,y,yhy,vec,
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      close(unit=299)

      
      
      
 1089   format(1x,A45)

 1090  format(1x, A15, 2x,I3)

 1091  format(1x,A5,2x,F12.8,A5,2x,F12.8,A5,2x,F12.8)
      
      
      
      
      
      
      end
