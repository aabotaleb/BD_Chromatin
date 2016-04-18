ccccccccccccccccccccccccccccccccccccccccc
ccc Modification: TONI (LHref)
c     Module that contains global parameters/variables

      module modglob
      
      implicit none

ccc   LHbeads
      ! nbNt:  beads in the N-term
      ! nbGh:  beads in the Glob-head
      ! nbCt:  beads in the C-term
      ! nbLH:  beads in the LH    
      integer nbNt,nbGh,nbCt,nbLH
      
ccc   LH coord and charges (arrays)
      ! rxyzNt,rxyzGh,rxyzCt: coordinates
      ! rqNt,rqGh,rqCt: charges
      double precision,dimension(:,:),allocatable :: rxyzNt
      double precision,dimension(:,:),allocatable :: rxyzGh
      double precision,dimension(:,:),allocatable :: rxyzCt
      double precision,dimension(:,:),allocatable :: rxyzLH
      double precision,dimension(:,:),allocatable :: rxyzLHeq
      double precision,dimension(:),allocatable :: rqNt,rqGh,rqCt
      double precision,dimension(:),allocatable :: rqLH
      

ccc   LH van der waals paremters
      ! revd_hh: histone-hitone (all LH beads)
      ! revd_hc: histone-core
      ! revd_hl: histone-linker
      ! revd_hhNt: histone-histone for Nt beads
      ! and so on
      double precision,dimension(:),allocatable :: revd_hh
      double precision,dimension(:),allocatable :: revd_hc
      double precision,dimension(:),allocatable :: revd_hl
      double precision,dimension(:),allocatable :: revd_ht
      double precision,dimension(:),allocatable :: revd_hhNt
      double precision,dimension(:),allocatable :: revd_hcNt
      double precision,dimension(:),allocatable :: revd_hlNt
      double precision,dimension(:),allocatable :: revd_htNt
      double precision,dimension(:),allocatable :: revd_hhGh
      double precision,dimension(:),allocatable :: revd_hcGh
      double precision,dimension(:),allocatable :: revd_hlGh
      double precision,dimension(:),allocatable :: revd_htGh
      double precision,dimension(:),allocatable :: revd_hhCt
      double precision,dimension(:),allocatable :: revd_hcCt
      double precision,dimension(:),allocatable :: revd_hlCt
      double precision,dimension(:),allocatable :: revd_htCt
      double precision,dimension(:),allocatable :: revd_hhLH
      double precision,dimension(:),allocatable :: revd_hcLH
      double precision,dimension(:),allocatable :: revd_hlLH
      double precision,dimension(:),allocatable :: revd_htLH

ccc   LH elastic parameters(arrays)
      ! rkstr*: streching constants
      ! rkben*: bending constants
      ! rktor*: torsion constants
      ! rkstreq: streching equilibrium (distance)
      ! rkbeneq: bending equilibrium (angle)
      ! rktoreq: torsion equilibrium (angle)
      ! lengthLH*: interbead LH distances (streching)
      ! betaLH*: interbead LH beta angle (bending)
      double precision,dimension(:),allocatable :: rkstr,rkben,rktor
      double precision,dimension(:),allocatable :: rkstrd2
      double precision,dimension(:),allocatable :: rkbend2
      double precision,dimension(:),allocatable :: rktord2
      double precision,dimension(:),allocatable :: rkstreq
      double precision,dimension(:),allocatable :: rkbeneq
      double precision,dimension(:),allocatable :: rktoreq
      double precision,dimension(:),allocatable :: lengthLH
      double precision,dimension(:),allocatable :: betaLH
      double precision,dimension(:),allocatable :: rkstrNt
      double precision,dimension(:),allocatable :: rkstrGh
      double precision,dimension(:),allocatable :: rkstrCt
      double precision,dimension(:),allocatable :: rkstrLH
      double precision,dimension(:),allocatable :: rkbenNt
      double precision,dimension(:),allocatable :: rkbenGh
      double precision,dimension(:),allocatable :: rkbenCt
      double precision,dimension(:),allocatable :: rkbenLH
      double precision,dimension(:),allocatable :: rktorNt
      double precision,dimension(:),allocatable :: rktorGh
      double precision,dimension(:),allocatable :: rktorCt
      double precision,dimension(:),allocatable :: rktorLH

ccc  LH connectivity (next bead in the sequence) and bending list
      integer,dimension(:),allocatable :: connLH
      integer,dimension(:),allocatable :: h_conn
      integer,dimension(:),allocatable :: bendLH
      integer,dimension(:),allocatable :: h_bend
      double precision,dimension(:),allocatable :: LHbound_randarray
      
ccc   LH elastic energy terms
      ! ELHstr: streching
      ! ELHben: bending
      ! ELHtor: torsion
      double precision ELHstr,ELHben,ELHtor

ccc   LH groups
      integer,dimension(:),allocatable :: LH_grp

ccc   Excluded volume parameters
      double precision k_exhl !!! Prefactor Ev LH/linkerDNA

ccc   Bound LH flag
      integer,dimension(:),allocatable :: LHbound
      integer,dimension(:),allocatable :: LHboundbds


      end module modglob
ccccccccccccccccccccccccccccccccccccccccc
