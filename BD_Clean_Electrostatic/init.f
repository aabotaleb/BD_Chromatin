C 	Main Subroutine that initializes coords, sets parameters, 
C	and implements MC loop

      subroutine init(nlb,n_c,nc3,n,n3,t_n,t_n3,Nq,Nq3, 
     + h_n,h_n3,ierr,myid,np,nbn)
              !Adding t_n3 parameter - Abotleb
      
      use modglob
      use mpi
      implicit NONE

c RCG: Definition of the dimension parameters
c     nbn     number of DNA beads
c     n_c     number of nucleosome cores
c     nc3     total number of nucleosome core position DoF
c     n       total number of beads in main chain: linker DNA beads + nucleosomes + addional flanking DNA beads
c     n3      total number of elements of the 3D (x,y,z) coors (nucleosome + linker DNA bead + flanking DNA) vectors 
c     t_n     total number of histone tail beads 
c     t_n3    total number of histone tail beads 3D coors -Abotaleb-
c     Nq      total number of pseudo charges per nucleosome core. 
c     Nq3     number of position DoF of these charges whithin a nucleosome core.
c     h_n     total number of linker histone beads 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


CCCCCCCCCCCCCCCCCCC   4/oct/2015 
c force(n3)	Translational gradient of the potential
c
c torque(n3)	Rotational gradient of the potential


C	1.Declare Variables!
      integer n_c, nc3, n, n3, t_n, steps,nlb,nbn(n_c)
      integer t_n3,h_n3   ! Abotaleb
      integer Nq, Nq3, freq, h_n,ih   
      integer ntpc, N_unit, n_att(5), n_succ(5), trialID     
      integer type(n), ki, m1, m2, m3
      integer t_grp(t_n), t_fix(t_n)           
      integer i,j,k,l, kk,ll, nm1, jm, j1,j2,j3, jm1, Nfort
      integer i1,i2,i3, N_flimit
      integer k_update
      integer startseed, seed
      integer ETIME
      integer facc
      integer laststep
      integer myid,ierr,np,restart,Ntrials,first,last
      integer Ngrids, occmatrix(14,14,10)
	  integer flpr
      integer flMg,fnonparLH
      integer firstn3,lastn3
      integer k1start,k1end
      integer k2start,k2end
      integer nCtNt ! number of N-term and C-term beads in one LH
      integer LH_nCtNt ! total number of N-term and C-term beads
      integer LHb
      integer k0,k1,k2,k3
      integer mattLH,maccLH,freqmLH,stepmLHlim
      integer signmLH,signmLHold
      integer LHcount
      integer LHnum
      integer nseq,nconn,nconn_old,nbend
      integer jmin,jmax,zz
      integer modeLHc          !!! Mode for the LH concentration algorithm



      
      character inputfile*50,corefile*50,tailfile*50,LH_input*50
      character LH_equil*50  
      character*10 Namefort
      character*10 oldfile

      
      REAL TIMEEND, TIMESTART, TIMEUSER, TIMESYSTEM
      REAL TARRAY(2), moveprob(4)

      common /random/ seed
      common /coeff/ debyell,fnonparLH !!! Common block for coefficients
      common /Ev/ Ev_LL,Ev_CC,Ev_TT,Ev_LC,Ev_CT,Ev_TL,Ev_TT1,Ev_TT2,
     +     Ev_CT1,Ev_CT2,Ev_HH,Ev_HT,Ev_HL1,Ev_HL2,Ev_HC
      logical	withlink,MOVEACCEPT,VERBOSE
      
      real*8, dimension(h_n)::drold,drnew,daold,danew,dbold,dbnew
      real*8, dimension(h_n)::dcold, dcnew
     
      double precision kbt, per, pi, dt
      double precision k_e
      double precision lo, ro, d1, go
      double precision debye
      double precision Cs, T
      double precision t_exv_e, t_exv_d
      double precision vdw_cut, Rcut, t_bcr
      double precision evd_tc, evd_tl, evd_cc, evd_cl, evd_ll
      double precision evd_hcG,evd_hlG,evd_hcC,evd_hlC, evd_hh
      double precision h_tc, Eb_tc, oldr, pacc      
      double precision gridsize,boundx,boundy,boundz,evd_link 

      double precision r(n3), a(n3), b(n3), c(n3)  !vectors for core-bead positions
      double precision aold(n3),bold(n3),cold(n3)
      double precision r_n(n3), a_n(n3), b_n(n3), c_n(n3)
      double precision a_dna(nc3),b_dna(nc3),c_dna(nc3)
      double precision a_dna_n(nc3),b_dna_n(nc3),c_dna_n(nc3)
      
      double precision alpha(n),beta(n),gamma(n) !stretching angles
      double precision alpha_n(n),beta_n(n),gamma_n(n)
      double precision alpha_p(n_c),beta_p(n_c),gamma_p(n_c)
      double precision alpha_p_n(n_c),beta_p_n(n_c),gamma_p_n(n_c)
      
      double precision length(n),length_n(n) !array holding distances
      double precision h, g, s, hd2, gd2, sd2, k_ex
      double precision phi_o(n_c), deltaphi, dist
      
      double precision t_chg(t_n), t_rad(t_n), t_mass(t_n) ! tail params
      double precision t_X(t_n), t_Y(t_n), t_Z(t_n)
      double precision t_X0(t_n), t_Y0(t_n), t_Z0(t_n)
      double precision t_X_n(t_n), t_Y_n(t_n), t_Z_n(t_n)
      double precision h_X(h_n),h_Y(h_n),h_Z(h_n), h_chg(h_n) !histone vars
      double precision h_X_n(h_n),h_Y_n(h_n),h_Z_n(h_n)
      double precision h_X0(3),h_Y0(3),h_Z0(3)
      double precision t_bond(t_n), t_bond_v(t_n), t_bond_c(t_n)  !tail value or constant
      double precision t_angle(t_n), t_angle_v(t_n), t_angle_c(t_n) 
      double precision t_Eb, t_Ea, t_Ec, t_Ev						! energy of tails
      
      double precision q_l, core_pos(Nq3),core_q(Nq) !arrays holding DiSCO charge coords
      double precision Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL !Energy contrib arrays
      double precision Ec_TT1, Ec_TT2, Ec_CT1, Ec_CT2  
      double precision Ec_HH, Ec_HT,Ec_TH,Ec_HL1, Ec_HL2, Ec_HC
      double precision E(6), delta1, delta2, delta3    
      double precision z(3), extra, randarray(17), qH1Gf, qH1Cf1,qH1Cf2 
      double precision temporary, Eold, Enew, Etot,EtotLH 
      double precision qH1G,qH1C1,qH1C2
     
      double precision Eoldpot,Enewpot,sum_Etot
      double precision dEpot,dEpar
      double precision rnew(n3),rold(n3)
      double precision mAc(3,3),mAcT(3,3),mAcI(3,3)
      double precision ac(3),bc(3),cc(3),ccc(3)
      double precision ac_old(3),bc_old(3),cc_old(3)
      double precision norm,proj
      double precision ampl,ranu
      
      double precision,allocatable :: yran(:),rvec(:),vec(:)
      
      double precision r_tmp(n3),a_tmp(n3), b_tmp(n3), c_tmp(n3)


      double precision debyell
      double precision Ev_LL, Ev_CC, Ev_TT, Ev_LC, Ev_CT, Ev_TL !Excluded volume 
      double precision Ev_TT1, Ev_TT2, Ev_CT1, Ev_CT2
      double precision Ev_HH, Ev_HT, Ev_HL1, Ev_HL2, Ev_HC
      double precision pot_send_array(25),sum_pot_send_array(25)
      double precision Esum(6)  

      integer e_counter,clock  
      double precision rx,ry,rz,rq

      double precision rehh,rehc,rehl,reht
      double precision kstr,kben,ktor
      double precision dx,dy,dz,dr,dr2
      double precision rxyz1(3),rxyz2(3),rxyz3(3)
      double precision rxyz12(3),rxyz23(3)
      double precision norm12,norm23,dotpr
      double precision cosbet,beta_tmp
      double precision delta

      double precision lengthLH_n(h_n),betaLH_n(h_n)
      double precision l1,l2
      double precision dEparpot
      double precision probLH,accprobLH ! probability for LH moves
      double precision deltaLH0,deltaLH ! typical displacement for a LH move
      double precision ratmLH0,ratmLH,ratmLHold,drat,dratold,ratmthr 
      double precision ddelLH,ddelLHth
      double precision h_X_old(h_n),h_Y_old(h_n),h_Z_old(h_n)
      double precision lengthLH_old(h_n),betaLH_old(h_n)
      double precision prek_exhl

      double precision LHconc  !!! Concentration of LHs/nucleosome
      double precision Pdice
      double precision randarray_regrow(2),twist_rand
 
      double precision force(n3), torque(n3)
      !Abotaleb added to account for forces isolated for each DNA bead isolated for each type
      double precision force_dnaStr_X(n),force_dnaStr_Y(n)
      double precision force_dnaStr_Z(n)
      double precision force_dnaBend_X(n),force_dnaBend_Y(n)
      double precision force_dnaBend_Z(n)
      double precision force_dnaTwist_X(n),force_dnaTwist_Y(n)
      double precision force_dnaTwist_Z(n)
      !Abotaleb added to isolate the energies 
      double precision ssle_Vec(n),ssbe_Vec(n),sstw_Vec(n)  ! Abotaleb Added to isolate the energy associated with each DNA Core/Bead
 
      ! Abotaleb  added to account for the forces due to tail beads
      double precision force_t(t_n3),force_tStrInt_X(t_n)
      double precision force_tStrInt_Y(t_n),force_tStrInt_Z(t_n)
      double precision force_tBend_X(t_n),force_tBend_Y(t_n)
      double precision force_tBend_Z(t_n)
      double precision t_Eb_num,t_Ea_num,Eb_tc_num ! numerical values of stretching and bending integrating the forces
     
      ! Abotaleb Electrostatic forces [April 2016 - Isolate the forces]
	  double precision force_c_LL(n3),force_c_LC(n3),force_c_CC(n3)
	  double precision force_c_HH(h_n3),force_c_HL1(h_n3),force_c_HL2(h_n3)
	  double precision force_c_HC(h_n3)
	  double precision force_c_TT1(t_n3),force_c_TT2(t_n3),force_c_TT(t_n3)
	  double precision force_c_TH(t_n3)
	  double precision force_c_TL(t_n3),force_c_CT1(t_n3),force_c_CT2(t_n3)
	  double precision force_c_CT(t_n3) 
      double precision force_c_TL_DNA(n3),force_c_CT1_DNA(n3)
	  double precision force_c_CT2_DNA(n3),force_c_CT_DNA(n3)
	  double precision force_c_TH_Histone(h_n3)
	  double precision force_c_HL1_DNA(n3),force_c_HL2_DNA(n3)
	  double precision force_c_HC_psdoChg_Core(Nq3)

      integer		   anneal_counter, decay_counter
      double precision anneal_temp, anneal_per, anneal_decay
      integer 		   anneal_flag
      logical gen_seed,gen_twist        

      
      ! Abotaleb Added to account for numerical diffentiation function calculated forces
      double precision force_dnaStr_X_num(n),force_dnaStr_Y_num(n)
      double precision force_dnaStr_Z_num(n)
      double precision force_dnaBend_X_num(n),force_dnaBend_Y_num(n)
      double precision force_dnaBend_Z_num(n)
      double precision force_dnaTwist_X_num(n),force_dnaTwist_Y_num(n)
      double precision force_dnaTwist_Z_num(n)

      parameter(k_e = 0.4151d0) ! Equilibrium DNA Bead Length
      parameter(pi = 3.14159265358979d0)
      parameter(lo=3.0d0)	!Debye-Huckel Params
      parameter(ro=4.8d0)	!Debye-Huckel Params
      parameter(d1=1.8d0)	!Debye-Huckel Params
      parameter(go=1.20d0*1.57079632679490d0) !Debye-Huckel Params
      parameter(t_exv_d = 1.8d0) !LJ radius tails
      parameter(t_exv_e = 0.1d0) !LJ prefactor 
      parameter(evd_tc = 1.8d0)	 !LJ tail-core
      parameter(evd_cc = 1.2d0)  !LJ core-core
      parameter(evd_cl = 2.4d0)  !LJ core-linker
      parameter(evd_ll = 3.6d0) !LJ linker-linker
      parameter(evd_tl = 2.7d0) !LJ tail-linker


C	2.Read Input!
C		if (myid.eq.0) then

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCC Debugging Purposes  Abdelrhman Modifications 8 Sept 2015         CCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c and  call GETARG(2, inputfile) replaced by   inputfile='input.run' 
c	  call GETARG(2, inputfile)
      inputfile='t1/input.run' 
      open (unit = 101, file = inputfile)
      write(*,'(A,i2)'),' INPUT FILE BEING READ BY.....',myid
      read(101,*)
      read(101,*) T, Cs, steps, restart, freq, (moveprob(i), i=1,4)
      read(101,*)
      read(101,*) t_bcr,delta1,delta2,delta3,Ntrials
      read(101,*)
      read(101,*) corefile,tailfile
      read(101,*)
      read(101,*) gen_seed,startseed    
      read(101,*)
      read(101,*) gen_twist,deltaphi
      read(101,*)
      read(101,*) flMg
      read(101,*)
      read(101,*) fnonparLH
      read(101,*)
      read(101,*) probLH,deltaLH
      read(101,*)
      read(101,*) prek_exhl
      read(101,*)
      read(101,*) LHconc,modeLHc,LHnum
      read(101,*)
      read(101,*) VERBOSE
      read(101,*) 
      read(101,*) anneal_flag,anneal_per,anneal_decay 
      close(101)
	  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       
       if (gen_seed) then
       call system_clock(COUNT=clock)
       startseed = clock
       end if

       call ranset(startseed)

       if (gen_twist.and.(.not.restart)) then
         if (myid.eq.0) then
         twist_rand = ranu()
         twist_rand = ranu()
C          print*,'random number for twist',twist_rand
          if (twist_rand.ge..6) then
          deltaphi=.3159
          else if (twist_rand.ge..3) then
          deltaphi=0.003
          else 
          deltaphi=-.3159
          end if
         end if
        end if

        call MPI_BCAST(deltaphi, 1, MPI_DOUBLE_PRECISION,
     +     0,MPI_COMM_WORLD,ierr)



! Automatic calculation of mean DNA twist according to DNA linker length
c The following phi_o are defined for lo=3nm (as in the PNAS) (from Rosana's vNRL')
       do i=1,n_c
       if(nbn(i).lt.2)then
       stop 'number of beads less than 2'
      elseif(nbn(i).eq.2)then
         phi_o(i)=0.9007+deltaphi/3
      elseif(nbn(i).eq.3)then
         phi_o(i)=-0.6754+deltaphi/4
      elseif(nbn(i).eq.4)then
         phi_o(i)=-0.3519+deltaphi/5
      elseif(nbn(i).eq.5)then
         phi_o(i)=-0.1463+deltaphi/6
      elseif((nbn(i).eq.6).or.(mod(nbn(i),6).eq.0))then
         phi_o(i)=(0.0215+deltaphi)/(nbn(i)+1)
      elseif(nbn(i).eq.7)then
         phi_o(i)=0.1178+deltaphi/8
      elseif(nbn(i).eq.8)then
         phi_o(i)=0.2025+deltaphi/9
      else
         stop 'number of beads greater than 8, and not multiple of 6'
      endif
      enddo


C		endif
C Initialize the main parameters for the LH MC move
      accprobLH = moveprob(1) + moveprob(2) + moveprob(3) + probLH
      deltaLH0 = deltaLH   
      freqmLH = 1000
      stepmLHlim = 1E6
      ddelLH = 0.05
      ddelLHth = 0.25
      ratmthr = 0.1
      ratmLH0 = 0.5
      ratmLHold = 0
      signmLHold = +1
      mattLH = 0
      maccLH = 0

C       withlink is now hardwired!
        withlink = T

C      if ((withlink).and.(myid.eq.0)) then
      if (withlink) then         !...............................................................WITHLINK LOOP START
         evd_hcG = 2.2d0        ! The evd_h* values are not consisten with previous papers,
         evd_hcC = 2.4d0        ! however, these are the values that had been used by Arya and Ogi.
         evd_hlG = 3.4d0        !
         evd_hlC = 3.6d0        !
         evd_hh  = 2.0d0		  ! Need to figure out what these are and find better values

         allocate(rxyzNt(nbNt,3))
         allocate(rxyzGh(nbGh,3))
         allocate(rxyzCt(nbCt,3))
         allocate(rqNt(nbNt))
         allocate(rqGh(nbGh))
         allocate(rqCt(nbCt))
         allocate(revd_hhNt(nbNt),revd_hcNt(nbNt),revd_hlNt(nbNt))
         allocate(revd_hhGh(nbGh),revd_hcGh(nbGh),revd_hlGh(nbGh))
         allocate(revd_hhCt(nbCt),revd_hcCt(nbCt),revd_hlCt(nbCt))
         allocate(revd_htNt(nbNt),revd_htGh(nbGh),revd_htCt(nbCt))
         allocate(revd_hh(h_n),revd_hc(h_n),revd_hl(h_n),revd_ht(h_n))
         allocate(rkstr(h_n),rkben(h_n),rktor(h_n))
         allocate(rkstrd2(h_n),rkbend2(h_n),rktord2(h_n))
         allocate(rkstreq(h_n),rkbeneq(h_n),rktoreq(h_n))
         allocate(lengthLH(h_n),betaLH(h_n))
         allocate(rkstrNt(nbNt),rkstrGh(nbGh),rkstrCt(nbCt))
         allocate(rkbenNt(nbNt),rkbenGh(nbGh),rkbenCt(nbCt))
         allocate(rktorNt(nbNt),rktorGh(nbGh),rktorCt(nbCt))
         allocate(LH_grp(h_n))
         allocate(rxyzLH(nbLH,3),rqLH(nbLH))
         allocate(rxyzLHeq(nbLH,3))
         allocate(revd_hhLH(nbLH),revd_hcLH(nbLH))
         allocate(revd_hlLH(nbLH),revd_htLH(nbLH))
         allocate(rkstrLH(nbLH),rkbenLH(nbLH),rktorLH(nbLH))
         allocate(connLH(nbLH))
         allocate(h_conn(h_n))
         allocate(bendLH(nbLH))
         allocate(h_bend(h_n))
         allocate(LHbound_randarray(h_n))

         nCtNt = (nbNt + nbCt) !Number of NTERM/CTERM LH Beads
         LH_nCtNt = n_c*nCtNt  !


ccc   Read 'LH.in'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCC Debugging Purposes  Abdelrhman Modifications 8 Sept 2015         CCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c and  call GETARG(3,LH_input)   replaced by   
c      call GETARG(3,LH_input)
      LH_input='t1/LH_N0G6C22.in'
         open(unit=101,file=LH_input)
         write(*,'(A,i2)'),' LH INPUT FILE BEING READ BY..',myid
         read(101,*)
ccc   N-term 
         read(101,*)
         do i = 1,nbNt
            read(101,*) rx,ry,rz,rq,rehh,rehc,rehl,reht,kstr,kben,ktor,
     +      nseq,nconn
            !!! Domain
            rxyzNt(i,:) = (/rx,ry,rz/)
            rqNt(i) = rq
            revd_hhNt(i) = rehh
            revd_hcNt(i) = rehc
            revd_hlNt(i) = rehl
            revd_htNt(i) = reht
            rkstrNt(i) = kstr
            rkbenNt(i) = kben
            rktorNt(i) = ktor

            !!! LH
            rxyzLH(nseq,:) = (/rx,ry,rz/)
            rqLH(nseq) = rq
            revd_hhLH(nseq) = rehh
            revd_hcLH(nseq) = rehc
            revd_hlLH(nseq) = rehl
            revd_htLH(nseq) = reht
            rkstrLH(nseq) = kstr
            rkbenLH(nseq) = kben
            rktorLH(nseq) = ktor
            connLH(nseq) = nconn
         enddo
C         print*,'FINISHED READING N-term'
ccc   G-head 
         read(101,*)
         do i = 1,nbGh
            read(101,*) rx,ry,rz,rq,rehh,rehc,rehl,reht,kstr,kben,ktor,
     +      nseq,nconn
            !!! Domain
            rxyzGh(i,:) = (/rx,ry,rz/)
            rqGh(i) = rq
            revd_hhGh(i) = rehh
            revd_hcGh(i) = rehc
            revd_hlGh(i) = rehl
            revd_htGh(i) = reht
            rkstrGh(i) = kstr
            rkbenGh(i) = kben
            rktorGh(i) = ktor

            !!! LH
            rxyzLH(nseq,:) = (/rx,ry,rz/)
            rqLH(nseq) = rq
            revd_hhLH(nseq) = rehh
            revd_hcLH(nseq) = rehc
            revd_hlLH(nseq) = rehl
            revd_htLH(nseq) = reht
            rkstrLH(nseq) = kstr
            rkbenLH(nseq) = kben
            rktorLH(nseq) = ktor
            connLH(nseq) = nconn
C                        print*,'nconn',nconn
           
         enddo
C         print*,'FINISHED READING G-head'

ccc   C-term 
         read(101,*)
         do i = 1,nbCt
            read(101,*) rx,ry,rz,rq,rehh,rehc,rehl,reht,kstr,kben,ktor,
     +      nseq,nconn
            !!! Domain
            rxyzCt(i,:) = (/rx,ry,rz/)
            rqCt(i) = rq
            revd_hhCt(i) = rehh
            revd_hcCt(i) = rehc
            revd_hlCt(i) = rehl
            revd_htCt(i) = reht
            rkstrCt(i) = kstr
            rkbenCt(i) = kben
            rktorCt(i) = ktor

            !!! LH
            rxyzLH(nseq,:) = (/rx,ry,rz/)
            rqLH(nseq) = rq
            revd_hhLH(nseq) = rehh
            revd_hcLH(nseq) = rehc
            revd_hlLH(nseq) = rehl
            revd_htLH(nseq) = reht
            rkstrLH(nseq) = kstr
            rkbenLH(nseq) = kben
            rktorLH(nseq) = ktor
            connLH(nseq) = nconn  
C                 print*,'nconn',nconn
         
         enddo
         close(101)

C         print*,'FINISHED READING C-term'
C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C         print*,'After MPI_BARRIER'

         !!! LH beads associated to a bending term
         !!! NOTE I haven't the foggiest idea why but the following lines sometimes cause a hang only on processor- something to do with MPI libraries
C         if (myid.eq.0) then


      
         !!! LH beads associated to a bending term
         nconn_old = 0
         do i=1,nbLH
            nconn = connLH(i)
            if((nconn==1).AND.(nconn_old==1)) then
               bendLH(i) = 1
            else
               bendLH(i) = 0
            endif
            nconn_old = nconn
C            write(*,*)'bendLH',bendLH(i),'connLH(i):',connLH(i)
         enddo
         
        do j=1,h_n
         do i=1,nbLH
           h_conn(j) = connLH(i)
           h_bend(j) = bendLH(i)
C           print*,h_conn(j),h_bend(j)
         enddo
        enddo
         
         
ccc Check LH
C		write(*,*)'connLH_dump',connLH

C      do i=1,nbLH
C            rx = rxyzLH(i,1)
C            ry = rxyzLH(i,2)
C            rz = rxyzLH(i,3)
C            rq = rqLH(i)
C
C            rehh = revd_hhLH(i)
C            rehc = revd_hcLH(i)
C            rehl = revd_hlLH(i)
C            reht = revd_htLH(i)
C            kstr = rkstrLH(i)
C            kben = rkbenLH(i)
C            ktor = rktorLH(i)
C            nconn = connLH(i)
C            nbend = bendLH(i)
C            !write(*,*) 'LH',i,nconn,nbend
C            write(*,*) 'LH',i,rx,ry,rz,rq,rehh,rehc,rehl,reht,
C     +           kstr,kben,ktor,nconn,nbend
C     		write(*,*)i,connLH(i)
C      enddo
      

         
C         call MPI_BCAST(bendLH,h_n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
ccc   Read 'LH_equil.in'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCC Debugging Purposes  Abdelrhman Modifications 8 Sept 2015         CCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
c and  call GETARG(4,LH_equil)   replaced by   LH_equil='LH_equil.in'
c		 call GETARG(4,LH_equil)
         LH_equil='t1/LH_N0G6C22_equil.in'
         open(unit=201,file=LH_equil)
         read(201,*)
         write(*,'(A,i2)'),' LH_EQUIL FILE BEING READ BY..',myid

ccc   N-term 
         read(201,*)
         do i = 1,nbNt
            read(201,*) rx,ry,nseq
            rkstreq(nseq) = rx
            rkbeneq(nseq) = ry
         enddo
ccc   G-head 
         read(201,*)
         do i = 1,nbGh
            read(201,*) rx,ry,nseq
            rkstreq(nseq) = rx
            rkbeneq(nseq) = ry
         enddo
ccc   C-term 
         read(201,*)
         do i = 1,nbCt
            read(201,*) rx,ry,nseq
            rkstreq(nseq) = rx
            rkbeneq(nseq) = ry
         enddo
      rkbeneq = rkbeneq*pi/180
         close(201)

      else     !..................................................WITHLINK ELSE

         qH1G = 0.0d0
         qH1C1 = 0.0d0
         qH1C2 = 0.0d0

         evd_hcG = 0.0d0
         evd_hcC = 0.0d0
         evd_hlG = 0.0d0
         evd_hlC = 0.0d0
         evd_hh  = 0.0d0       
      endif   !.........................................................WITHLINK ENDIF

c	  write(*,*) "LH Parameters"	
c      write(*,*) 'Charges: Globular Head, C1, C2: ',qH1G,qH1C1,qH1C2
c      write(*,*) 'LJ:,',evd_hcG,evd_hcC,evd_hlG,evd_hlC,evd_hh
c      write(*,*) 'Nterm f_c str, ben, tor: ',rkstrNt,rkbenNt,rktorNt
c      write(*,*) 'GH f_c str, ben, tor:',rkstrGh,rkbenGh,rktorGh
c      write(*,*) 'CTER force constants str, ben, tor:',rkstrCt,rkbenCt,rktorCt

C	print*,'Finished with LH_equil!',myid

      qH1Gf = qH1G
      qH1Cf1 = qH1C1
      qH1Cf2 = qH1C2
      
!parameters
      h=100.d0
      h_tc=100.d0
      N_unit = n / n_c
      ntpc = t_n / n_c
      debye = 0.736d0 * dsqrt( (Cs/0.050d0) * (298/T) )

c     Parameters that depend on the presence of Mg
c     - DNA persistence length
c     - Debye length of linker DNA /linker DNA interactions
      if(flMg) then
	     per = 30.0d0
         debyell = 2.5
      else
         per = 50.0d0
         debyell = debye
      endif
cccccccccccc

      q_l = lo * (-5.8824d0) *
     +           (   1.901912136835033d-08 * (Cs*1000)**3 
     +             - 8.211102025728157d-06 * (Cs*1000)**2
     +             + 7.554037628581672d-03 * (Cs*1000)
     +             + 3.524292543853884d-01         )


      kbt = 0.5962d0 * T / 300.0d0
cccccccccccc
      k_ex = 0.001d0*kbt      

      k_exhl = k_ex*prek_exhl !!! prefactor Ev LH/linkerDNA

      h = h* kbt / lo / lo  
      hd2 = 0.5d0 * h
      g = per * kbt / lo
      gd2 = 0.5d0 * g
      s =  43.17d0 / lo
      sd2 = 0.5d0 * s
      h_tc = h_tc *kbt /lo /lo  
      Rcut = 4.0d0 
      
      do i = 1, 100
          temporary = k_e*q_l*q_l*dexp(- debye*Rcut ) / Rcut
          if ( temporary .lt. 5.0d-3 ) then
              go to 500
          else    
              Rcut = Rcut + 1.0d0
          end if
      end do
  500 vdw_cut = 4.0d0
      do i = 1, 10
          temporary = k_ex*((evd_ll/vdw_cut)**12-(evd_ll/vdw_cut)**6)
          if ( dabs(temporary) .lt. 5.0d-3 ) then
              go to 600
          else    
              vdw_cut = vdw_cut + 1.0d0
          end if
      end do
 
C 3. Write main simulation parameters to fort.13 and screen 
 	 
  600 continue 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	  If (myid.eq.0) then
         if (flMg) then
      write(*,*)'............................
     +................. SIMULATION WILL HAVE Mg2+'
         endif
      	if (VERBOSE) then
      write(*,*)'.................
     +...............................WRITING VERBOSE OUTPUT'
	    endif
	    if (anneal_flag.ne.0) then
	    write(*,*)'.................
     +.....................WILL PERFORM SIMULATED ANNEALING'
        end if
      End if 
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (myid.eq.0) then !....................................MP WRITE PARAMS to SCREEN AND FILE START
	  	      
      write(*, fmt=9007) 'n_c = ', n_c  
      write(*, fmt=9007) 'n = ', n  
      write(*, fmt=9007) 'N_unit = ', N_unit 
      write(*, fmt=9007) 'nptc = ', ntpc
      write(*, fmt=9006) 'T = ', T, 'K'
      write(*, fmt=9006) 'Cs = ', Cs, 'M'      
      write(*, fmt=9006) 'debye = ', debye, 'nm^-1'
      write(*, fmt=9006) 'q_l = ', q_l, 'e'
      write(*, fmt=9006) 't_bcr = ', t_bcr, ' '        
      write(*, fmt=9006) 'kbt = ', kbt, 'kcal/mol'      
      write(*, fmt=9006) 'k_ex = ', k_ex, 'kcal/mol'       
      write(*, fmt=9006) 'h = ', h, 'kcal/mol/nm^2'       
      write(*, fmt=9006) 'g = ', g, 'kcal/mol'       
      write(*, fmt=9006) 's = ', s, 'kcal/mol'       
      write(*, fmt=9006) 'h_tc = ', h_tc, 'kcal/mol/nm^2'     
      write(*, fmt=9006) 'delta_phi = ', deltaphi, 'rad'       
      write(*, fmt=9006) 'Rcut = ', Rcut, 'nm'      
      write(*, fmt=9006) 'vdw_cut = ', vdw_cut, 'nm' 
      write(*, fmt=9006) 't_exv_d = ', t_exv_d, 'nm'   
      write(*, fmt=9006) 't_exv_e = ',t_exv_e, 	'(kcal/mol)'
      write(*, fmt=9006) 'evd_tl = ', evd_tl, 'nm' 
      write(*, fmt=9006) 'evd_tc = ', evd_tc, 'nm'       
      write(*, fmt=9006) 'evd_cc = ', evd_cc, 'nm' 
      write(*, fmt=9006) 'evd_cl = ', evd_cl, 'nm' 
      write(*, fmt=9006) 'evd_ll = ', evd_ll, 'nm'
C      write(*, fmt=9006) 'evd_lnk= ', evd_link, 'nm'
      write(*, fmt=9006) 'evd_hcG = ', evd_hcG, 'nm'
      write(*, fmt=9006) 'evd_hcC = ', evd_hcC, 'nm'
      write(*, fmt=9006) 'evd_hlG = ', evd_hlG, 'nm'
      write(*, fmt=9006) 'evd_hlC = ', evd_hlC, 'nm'
      write(*, fmt=9007) 'steps = ', steps
      write(*, fmt=9007) 'restart = ',  restart       
      write(*, fmt=9007) 'freq = ', freq
      write(*, fmt=9010) 'corefile = ', corefile 
      write(*, fmt=9010) 'tailfile = ', tailfile
      write(*, fmt=9007) 'iseed = ', startseed
C      write(*, fmt=9012) 'withlink = ', withlink
C      write(*, fmt=9006) 'q_h1G = ', qH1G, 'e'
C      write(*, fmt=9006) 'q_h1C1 = ', qH1C1, 'e'         
C      write(*, fmt=9006) 'q_h1C2 = ', qH1C2, 'e'
              

	      
	  
	  write(unit=13, fmt=9007) 'n_c = ', n_c  
      write(unit=13, fmt=9007) 'N_unit = ', N_unit 
      write(unit=13, fmt=9007) 'nptc = ', ntpc
      write(unit=13, fmt=9006) 'T = ', T, 'K'
      write(unit=13, fmt=9006) 'Cs = ', Cs, 'M'      
      write(unit=13, fmt=9006) 'debye = ', debye, 'nm^-1'
      write(unit=13, fmt=9006) 'q_l = ', q_l, 'e'
      write(unit=13, fmt=9006) 't_bcr = ', t_bcr, ' '        
      write(unit=13, fmt=9006) 'kbt = ', kbt, 'kcal/mol'      
      write(unit=13, fmt=9006) 'k_ex = ', k_ex, 'kcal/mol'       
      write(unit=13, fmt=9006) 'h = ', h, 'kcal/mol/nm^2'       
      write(unit=13, fmt=9006) 'g = ', g, 'kcal/mol'       
      write(unit=13, fmt=9006) 's = ', s, 'kcal/mol'       
      write(unit=13, fmt=9006) 'h_tc = ', h_tc, 'kcal/mol/nm^2'     
      write(unit=13, fmt=9006) 'delta_phi = ', deltaphi, 'rad'       
      write(unit=13, fmt=9006) 'Rcut = ', Rcut, 'nm'      
      write(unit=13, fmt=9006) 'vdw_cut = ', vdw_cut, 'nm' 
      write(unit=13, fmt=9006) 't_exv_d = ', t_exv_d, 'nm'   
      write(unit=13, fmt=9006) 't_exv_e = ',t_exv_e, 
     +				       '(kcal/mol)'
      write(unit=13, fmt=9006) 'evd_tl = ', evd_tl, 'nm' 
      write(unit=13, fmt=9006) 'evd_tc = ', evd_tc, 'nm'       
      write(unit=13, fmt=9006) 'evd_cc = ', evd_cc, 'nm' 
      write(unit=13, fmt=9006) 'evd_cl = ', evd_cl, 'nm' 
      write(unit=13, fmt=9006) 'evd_ll = ', evd_ll, 'nm'
C      write(unit=13, fmt=9006) 'evd_lnk= ', evd_link, 'nm'
      write(unit=13, fmt=9006) 'evd_hcG = ', evd_hcG, 'nm'
      write(unit=13, fmt=9006) 'evd_hcC = ', evd_hcC, 'nm'
      write(unit=13, fmt=9006) 'evd_hlG = ', evd_hlG, 'nm'
      write(unit=13, fmt=9006) 'evd_hlC = ', evd_hlC, 'nm'
      write(unit=13, fmt=9007) 'steps = ', steps
      write(unit=13, fmt=9007) 'restart = ',  restart       
      write(unit=13, fmt=9007) 'freq = ', freq
      write(unit=13, fmt=9010) 'corefile = ', corefile 
      write(unit=13, fmt=9010) 'tailfile = ', tailfile
      write(unit=13, fmt=9007) 'iseed = ', startseed
C      write(unit=13, fmt=9012) 'withlink = ', withlink
C      write(unit=13, fmt=9006) 'q_h1G = ', qH1G, 'e'
C      write(unit=13, fmt=9006) 'q_h1C1 = ', qH1C1, 'e'         
C      write(unit=13, fmt=9006) 'q_h1C2 = ', qH1C2, 'e'
      close(unit=13)

	  endif	!.....................................................................MP WRITE PARAMS END
	  
C	  call MPI_Barrier(MPI_COMM_WORLD,ierr)
C	  write(*,'(A,i2,A)'),' HELLO FROM PROCESS ',myid,'!'
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	  
      nm1 = n-1
      
	  call MPI_Barrier(MPI_COMM_WORLD,ierr)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/201
C       1-    Read coordinate file          /////////////////////////////
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/2015

        call startconf(n_c,nbn,n,n3,type,r,ro,a,b,c,lo,np,myid,ierr)

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/201
C       2-    Read core file          /////////////////////////////
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/2015

      call readcore(corefile,Nq,Nq3,core_pos,core_q,myid)
      
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/201
C       3-    Read tail file          /////////////////////////////
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/2015

      
      call readtail(tailfile,n_c,n,t_n,t_grp,t_fix,t_chg,t_rad,t_mass,
     +     t_X,t_Y,t_Z, t_X0,t_Y0,t_Z0, t_bond,t_bond_v,t_bond_c, 
     +     t_angle, t_angle_v, t_angle_c,n3,r,a,b,c,h_n,h_X,h_Y,
     +     h_Z,h_chg,ro, qH1G, qH1C1, qH1C2,h_X0,h_Y0,h_Z0,withlink,
     +	   myid,nbn,type)
      
      do i = 1, t_n
          t_chg(i) = t_chg(i)*t_bcr
      end do
      
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/201
C       4-    Read grids file          /////////////////////////////
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/2015
      

C      if (myid.eq.0) then
! call grids(Ngrids,gridsize,boundx,boundy,boundz,occmatrix,myid)
C      end if 
      
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/201
C       5-    Write LH Equilbrium bond length and equilbrium bend angle/////////////////////////////
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 4/oct/2015


	  if ((myid.eq.0).and.withlink) then
      open(unit=101,file='LH.out',status='unknown')
      write(101,*) '# Bead      Equil bond [nm]      Equil bend [rad]'
      write(*,*),'LH.OUT SUCCESSFULLY WRITTEN'
      
      do i = 1,nbLH
         rx = rkstreq(i)
         ry = rkbeneq(i)
         write(101,*) i,rx,ry
      enddo
      close(101)
      end if 

      laststep = 0
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C      print*, 'starting update_mod for the 1st time'
        aold=a
        bold=b
        cold=c
         call update_mod(n_c,nc3,n,n3,type,r,ro,d1,go,a,b,c,
     +        alpha,beta,gamma,length,a_dna,b_dna,c_dna,
     +        alpha_p,beta_p,gamma_p,1)
     
C  The update_mod loops update a,b,c differently on diff processors so we synchronize them here
      call MPI_BCAST(r, n3, MPI_DOUBLE_PRECISION, 
     +	    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(a, n3, MPI_DOUBLE_PRECISION, 
     +	    0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(b, n3, MPI_DOUBLE_PRECISION, 
     +	    0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(c, n3, MPI_DOUBLE_PRECISION, 
     +	    0,MPI_COMM_WORLD,ierr)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
ccc Initial streching distances and bending angles
      if (withlink) then ! ......................................WITHLINK IF START
       do i=1,h_n
         !!! Stretching (bond length)
         lengthLH(i) = 0.
         nconn = h_conn(i)
         if(nconn==1) then
            k1 = i
            k2 = k1 + 1
            rxyz1 = (/h_X(k1),h_Y(k1),h_Z(k1)/)
            rxyz2 = (/h_X(k2),h_Y(k2),h_Z(k2)/)
            !!! Distance
            rxyz12 = rxyz2 - rxyz1
            dr2 = dot_product(rxyz12,rxyz12)
            dr = sqrt(dr2)
            lengthLH(k1) = dr
         endif
c         write(*,*) 'l_LH',i,lengthLH(i)
         !!! Bending (bond angle)
         betaLH(i) = 0.
         nbend = h_bend(i)
         if(nbend==1) then
            k1 = i - 1
            k2 = i
            k3 = i + 1
            rxyz1 = (/h_X(k1),h_Y(k1),h_Z(k1)/)
            rxyz2 = (/h_X(k2),h_Y(k2),h_Z(k2)/)
            rxyz3 = (/h_X(k3),h_Y(k3),h_Z(k3)/)
            !!! Angle
            rxyz12 = rxyz2 - rxyz1
            rxyz23 = rxyz3 - rxyz2
            dr2 = dot_product(rxyz12,rxyz12)
            norm12 = sqrt(dr2)
            dr2 = dot_product(rxyz23,rxyz23)
            norm23 = sqrt(dr2)
            dotpr = dot_product(rxyz12,rxyz23)
            cosbet = dotpr/(norm12*norm23)
            if(cosbet.gt.1) then
               beta_tmp = 0
            else if(cosbet.lt.-1) then
               beta_tmp = pi
            else
               beta_tmp = acos(cosbet)
            endif
            betaLH(k2) = beta_tmp
c            write(*,*) 'beta tmp',beta_tmp,cosbet,dotpr
         endif
c         write(*,*) 'b_LH',i,betaLH(i)         
       enddo
		
		
		!!Get some random numbers for LHboundsub
		if (myid.eq.0) then
		  do i = 1,h_n
		    LHbound_randarray(i) = ranu()
		  enddo
		endif 
		

	 	call MPI_BCAST(LHbound_randarray, h_n, MPI_DOUBLE_PRECISION, 
     +	    0,MPI_COMM_WORLD,ierr)

		
      call LHboundsub(myid,np,n_c,h_n,startseed,LHconc,modeLHc,LHnum,
     +ierr,VERBOSE)
     
      end if     !....................................................WITHLINK ENDIF

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       call  potential(n_c,n,n3,type,r,a,b,c,alpha,beta,gamma,
     +	  length,beta_p,lo,hd2,gd2,sd2,k_e,k_ex,Nq,Nq3,core_pos,core_q,
     +    debye,q_l,phi_o,E,t_n,t_chg,t_X, t_Y, t_Z,t_X0, t_Y0, t_Z0,
     +    h_X,h_Y,h_Z,h_chg,h_n,evd_hcG,evd_hcC,evd_hlG,evd_hlC,
     +    t_bond, t_bond_v, t_bond_c, t_angle, t_angle_v, t_angle_c, 
     +    t_exv_e, t_exv_d, t_Eb, t_Ea, t_Ec, t_Ev,t_grp, t_fix,
     +    vdw_cut, Rcut, evd_tc, evd_tl, evd_cc, evd_cl, evd_ll,
     +    Ec_LL, Ec_CC, Ec_TT, Ec_LC, Ec_CT, Ec_TL,Ec_TT1, 
     +    Ec_TT2, Ec_CT1, Ec_CT2, Ec_HH, Ec_HT, Ec_HC, Ec_HL1,
     +    Ec_HL2, h_tc, Eb_tc, withlink, evd_link, evd_hh,
     +    np,myid,ierr,1)
cccccccccccccccccccccccccccccccccc
        call gfat(n_c,nc3,n,n3,type,r,a,b,c,E,
     +                alpha,beta,gamma,length,
     +                a_dna,b_dna,c_dna,
     +                alpha_p,beta_p,gamma_p,
     +                lo,ro,d1,go,h,g,s,k_e,k_ex,hd2,gd2,sd2,
     +                Nq,Nq3,core_pos,core_q,
     +                debye,q_l,phi_o,force,torque,t_n,t_n3, ! t_n,t_n3 number of tails and coors
     +                force_dnaStr_X,force_dnaStr_Y,force_dnaStr_Z,
     +                force_dnaBend_X,force_dnaBend_Y,force_dnaBend_Z,
     +             force_dnaTwist_X,force_dnaTwist_Y,force_dnaTwist_Z,
     +             force_c_LL,force_c_LC,force_c_CC,             ! forces affect DNA Beads  
     +             force_c_HH,force_c_HL1,force_c_HL2,force_c_HC,               ! Linker Histone-Linker Histone interactions  !linker Histone - linker DNA [Adjacent- Non adjacent] Histone Side
     +  force_c_HL1_DNA,force_c_HL2_DNA,force_c_HC_psdoChg_Core,          ! Linker Histone-Linker Histone interactions  !linker Histone - linker DNA [Adjacent- Non adjacent] DNA Linker-Core Side
     +             force_c_TT1,force_c_TT2,force_c_TT,	 ! Tail internal forces - Tail Linker Histone interactions
     +              force_c_TH,force_c_TH_Histone,       ! Tail -Histone interactions [Tail Side - Histone Side]
     +             force_c_TL,force_c_CT1,force_c_CT2,force_c_CT,               ! Tail external forces (Tail linker - Tail Core)  Tail Side
     + force_c_TL_DNA,force_c_CT1_DNA,force_c_CT2_DNA,force_c_CT_DNA,           ! Tail external forces (Tail linker - Tail Core)  DNA Side
     +             Ec_LL ,Ec_LC ,Ec_CC,Ec_HH,Ec_HL1,Ec_HL2,Ec_HC,    ! Electrostatics Energies
     +             Ec_TT1,Ec_TT2,Ec_TT,Ec_TL,Ec_CT1,Ec_CT2,Ec_CT,    ! cont. Electrostatic energies                      
     +            h_X,h_Y,h_Z,h_chg,h_n,h_n3,
     +                t_grp,t_fix,t_chg,t_rad,t_mass,        ! tail business
     +                t_X,t_Y,t_Z, t_X0,t_Y0,t_Z0,           ! >>
     +                t_bond,t_bond_v,t_bond_c,              ! >>
     +                t_angle, t_angle_v, t_angle_c,force_t, ! >>
     +                force_tStrInt_X,force_tStrInt_Y,force_tStrInt_Z,
     +                force_tBend_X,force_tBend_Y,force_tBend_Z,    ! isolating internal stretching and bending forces
     +                 np,myid,ierr)               

      write(*, fmt=9006) 'Exact Stretch Energy in tail beads = ',
     +                          t_Eb, 'kCal'      
      write(*, fmt=9006) 'Exact bedning Energy in tail beads = ',
     +                          t_Ea, 'kCal'
      write(*, fmt=9006) 'Approx Stretch Energy in tail beads = ',
     +                          t_Eb_num, 'kCal'
      write(*, fmt=9006) 'Approx bending Energy in tail beads = ',
     +                          t_Ea_num, 'kCal'
      
        ! Adding t_n,t_n3,force_t for calculating forces associated with tail beads  - Abotaleb -
      
      
 9000 format(2x, E14.7, 2x, E14.7, 2x, E14.7)
 9001 format(I14)


 9002 format(2x, E14.7, 2x, E14.7, 2x, E14.7, 
     +       2x, E14.7, 2x, E14.7, 2x, E14.7,
     +       2x, E14.7, 2x, E14.7, 2x, E14.7, 2x, E14.7,
     +       2x, E14.7,
     +       2x, E14.7, 2x, E14.7, 2x, E14.7)

 9003 format(2x, E14.7, 2x, E14.7, 2x, E14.7, 
     +       2x, E14.7, 2x, E14.7, 2x, E14.7,
     +       4x, E14.7, 2x, E14.7, 4x, E14.7, 2x, E14.7)
 9004 format(2x, I14, 2x, E14.7)
 9005 format(2x, E14.7, 2x, E14.7, 2x, E14.7, 2x, E14.7)
 9006 format(A20, 2X, E14.7, 2X, A)
 9007 format(A20, 2X, I14)
 9008 format(2x, E14.7, 2x, E14.7)
 9009 format(2x, A18, 2x, I9,I4, 2x, A14)
 9010 format(A20,2X,A)
 9011 format(2x,A,2x,I2,2X,A,2X,I6,2X,A,2x,I5)
 9012 format(A20, 2X, L14)
 9014 format(2x, E14.7, 2x, E14.7, 2x, E14.7, 
     +       2x, E14.7, 2x, E14.7, 2x, E14.7,
     +       4x, E14.7, 2x, E14.7, 4x, E14.7, 2x, E14.7,
     +       2x, E14.7, 2x, E14.7, 2x, E14.7, 2x, E14.7, 2x, E14.7)
9020  format(E14.7)
      
 9999 end
