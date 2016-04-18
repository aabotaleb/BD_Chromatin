c***************************************************************
      subroutine testgh(n,xc,fc,gc,y,yhy,vec,
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
     +    tested_comp,analytic,
     +    np,myid,ierr)
c
c GOAL:    Test user-supplied gradient (G) and Hessian (H)
c ----     routines corresponding to a given function F.
c          Testing H is optional.
c
c METHOD:    Derivatives are tested using a Taylor expansion of F
c -------    around a given point XC. The Taylor series is expanded 
c      at XC + eps*Y where Y is a random  perturbation vector
c      and eps is a scalar. If we denote the dot product of 2
c      vectors A and B as (A,B), we can write our expansion as
c
c      F(XC+eps*Y) = F(XC) + eps * (G,Y) + 1/2*(eps**2) * (Y,HY)
c                          + O(eps**3),
c
c      where G and H are both evaluated at XC, and HY denotes a
c      Hessian/vector product. If only G routines are tested, the
c      second-order Taylor term is zero, and the truncation error
c      is O(eps**2).
c
c      Our test is performed by computing this Taylor approx. at
c      smaller and smaller values of eps and checking to see
c      whether correct truncation errors are obtained --
c      O(eps**2) and  O(eps**3) if the approx. is correct upto
c      the G and H terms, respectively.
c
c      We divide eps by 2 at every step and test if indeed our
c      truncation errors decrease in the rate they should.
c      (i.e., if the error corresponding to eps is E1,
c      the error for eps/2 should be E1/4 if the gradient
c      is correct, and E1/8 if the Hessian is also correct).
c      Our value "RATIO" computes this factor of the old/new
c      errors.
c
c OUTPUT:    A series of values for RATIO is printed for each eps
c -------    until the truncation error and/or eps is very small.
c      If RATIO tends to 4 or 8 as eps is decreased,
c      G  is correct or G&H are correct, respectively.
c      If RATIO tends to 2, which is O(eps), neither G nor
c      H are correct. (If the error is larger than O(eps),
c      it is likely that this testing routines has a bug ...)
c
c      Keep in mind that the reliable values of RATIO should
c      occur when: (1) eps is not too large and not too small,
c      and (2) the difference between the  F(XC+eps*Y) and the
c      Taylor series approximation is of reasonable magnitude.
c      (The values of eps and the errors appear in the output).
c      In other words, a very accurate value of RATIO should
c      appear around the middle of our series. If the user has
c      any doubts, a different starting point and/or
c      perturbation vector should be tried.
c
c USAGE:     First, the user must supply the following input
c ------     variables in the function call:
c
c      N      - dimension (number of variables for F)
c      XC(N)  - our current vector
c      FC     - the function value at XC
c      GC(N)  - the gradient vector at XC, on input.
c               (NOTE: GC may be changed - see the note above).
c      Y(N)   - a random perturbation vector (Y should be chosen
c               so that F(XC+Y) is in a reasonable range for the
c               problem)
c      YHY    - the matrix inner product -- (Y,HY) -- representing
c               the dot product of Y with the Hessian/vector
c               product, HY, where H is evaluated at XC.
c               (if only the gradient is tested, set YHY to zero).
c      VEC(N) - a work vector
c
c      Second, the user must modify this routine by replacing
c      the sample function call given here with the appropriate
c      call for his/her problem. The inserted routine call
c      (just before the '30 CONTINUE' statement) should
c      produce a new function value, FVEC, for each new vector
c      VEC=XC+eps*Y. For example, if the user's subroutine
c      OBJFCT(N,X,F,G) computes F and G at X, the
c      required insertion is: CALL OBJFCT(N,VEC,FVEC,GC).
c      (NOTE: The new gradient at VEC is not required, but
c      our GC vector can be used here).
c
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      integer n
      double precision fc, yhy, xc(n),gc(n),y(n),vec(n)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c inserted to use by update_mod and potential functions
c        Variables needed for the update_mod routine
      integer n_c,nc3,num,n3, type(n),myid,ierr,np,index,j1,j2,j3,h_n
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
      !     local variables needed by the function 
      integer i,i1,i2,i3
      double precision rs(n3),as(n3),bs(n3),cs(n3)
      double precision alphas(n),betas(n),gammas(n),lengths(n)
      double precision a_dnas(nc3),b_dnas(nc3),c_dnas(nc3)
      double precision alpha_ps(n_c),beta_ps(n_c),gamma_ps(n_c)
      double precision Es(6)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer tested_comp,analytic
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      character(len=1024) :: filename

      if (analytic.eq.1) then
          if (tested_comp.eq.1) then
          write (filename, "(A40)") 'TestRes_Stretching_Analytic.txt'
          else 
              if (tested_comp.eq.2) then
          write (filename, "(A40)") 'TestRes_Bending_Analytic.txt'
              else
                  if (tested_comp.eq.3) then
           write (filename, "(A40)") 'TestRes_Twisting_Analytic.txt'
                  endif
              endif
          endif
      else
          if (tested_comp.eq.1) then
            write (filename, "(A40)") 'TestRes_Stretching_Numerical.txt'
          else 
              if (tested_comp.eq.2) then
              write (filename, "(A40)") 'TestRes_Bending_Numerical.txt'
              else
                  if (tested_comp.eq.3) then
              write (filename, "(A40)") 'TestRes_Twisting_Numerical.txt'
                  endif
              endif
          endif       
      endif
      

       open(unit=1361,name=filename,access='SEQUENTIAL',
     +        status='unknown')
       
      one    = 1.d0
      half   = 0.5d0
      call mcheps(epsmch)
      epslim = epsmch * 1.d+1
      epsmin = epsmch * 1.d+10
c     epsmin = epsmch * 1.d+5
      eps    = one  
      mp     = 6

      write(unit=1361,fmt=901)
      gy = 0.d0
      do 1 i = 1,n
         gy = gy + gc(i)*y(i)
    1 continue
      
       write(unit=1361,fmt=902) fc,gy, yhy,epsmch
       write(unit=1361,fmt= 910)

   10 temp=diff

      do 15 i = 1,n
         vec(i) = xc(i) + eps*y(i)
   15 continue

      nout = 0
cccccccccccccccccccccccccccccccccc
cccc     our objective function is inserted here
      as=aold
      bs=bold
      cs=cold

      call update_mod(n_c,nc3,num,n3,type,vec,ro,d1,go,as,bs,cs,
     +        alphas,betas,gammas,lengths,a_dnas,b_dnas,c_dnas,
     +        alpha_ps,beta_ps,gamma_ps,100)
     
      
      call potential(n_c,num,n3,type,vec,as,bs,cs,alphas,betas,gammas,
     +	  lengths,beta_ps,lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,Es,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_TH, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    np,myid,ierr,100)
      fvec=Es(tested_comp)
cccccccccccccccccccccccccccccccccccccccccccccccc       
   30 continue
      
      taylor = fc + (eps*gy) + ( (eps**2) * half  * yhy )
      diff   = fvec - taylor

      if ( abs(diff).lt.(epslim) ) then
         write(unit=1361,fmt=904) epslim
         goto 50 
      endif

      ratio = temp / diff
      if (eps.eq.one) then
           write(unit=1361,fmt=912) eps,fvec,taylor,diff
      else
           write(unit=1361,fmt= 911) eps,fvec,taylor,diff,ratio
      endif

      eps = eps * half
      if (eps .lt. epsmin) goto 50

      goto 10

   50 return
  901 format(///t10,'ENTERING TESTGH ROUTINE:'///) 
  902 format(t5,  'The function value at X               = ',
     + 1PE16.8/t5,'The first-order Taylor term,  (G, Y)  = ',
     + 1PE16.8/t5,'The second-order Taylor term, (Y,HY)  = ',
     + 1PE16.8//t5,'The computed machine precision        = ',
     + 1PE16.8//) 
  904 format(/t5,'DIFF is very small (LESS THAN ', 1PE16.8,')'/)
  910 format(4x,'EPS',10x,' F   ',10x,' TAYLOR',9x,
     +  ' DIFF.',12x,'RATIO'/)
  911 format(1x,F14.7,1x,F14.7,1x,F14.7,1x,F14.7,1x,F14.7)
  912 format(1x,F14.7,F14.7,F14.7,F14.7)
      close(unit=1361)
      end
C************************************************************
      subroutine mcheps(eps)
      double precision eps,one,two,half
      one  = 1.d0
      two  = 2.d0
      half = one / two
      eps  = one
    1 eps  = eps * half
      if ( (one + eps) .eq. one) goto 2
      goto 1
    2 continue
      eps  = two * eps
      return 
      end
