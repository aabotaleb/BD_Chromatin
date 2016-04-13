!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Modification: TONI (LHref)
!!! This subroutine updates the elastic distances and angles associated
!!! to the selected bead LHb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine updateLH(h_X,h_Y,h_Z,length,beta,&
     h_n,LHb,pi)


  !!! Variables
  use modglob !!! We activate the global parameters/variables

  implicit NONE

  integer h_n
  double precision length(h_n),beta(h_n)
  double precision h_X(h_n),h_Y(h_n),h_Z(h_n)
  integer LHb
  double precision xd(3)
  double precision rold(3),rnew(3)
  integer i,ic,i0,i1,i2,ilow,ihigh
  integer j1,j2,j3
  integer k1,k2,k3
  double precision rxyz1(3),rxyz2(3),rxyz3(3),drxyz(3)
  double precision rxyz12(3),rxyz23(3)
  double precision dr2
  double precision norm12,norm23
  double precision dotpr,cosbet,beta_tmp
  double precision pi
  integer nconn,nbend
  integer imin,imax

  !!! Stretching (bond length): It can affect two bonds
  !!! Previous and next bonds
  imin = LHb - 1
  imax = LHb
  do i=imin,imax
     nconn = h_conn(i)
     if(nconn==1) then
        j1 = i
        j2 = i + 1
        rxyz1(:) = (/h_X(j1),h_Y(j1),h_Z(j1)/)
        rxyz2(:) = (/h_X(j2),h_Y(j2),h_Z(j2)/)
        drxyz = rxyz2 - rxyz1
        dr2 = dot_product(drxyz,drxyz)
        length(i) = sqrt(dr2)
        !write(*,*) 'lgth',j1,j2,length(i)
     endif
  enddo

  !!! Bending (bond angle): It can affect three bonds
  !!! Previous/Current/Next bead angles
  imin = LHb - 1
  imax = LHb + 1
  do i=imin,imax
     nbend = h_bend(i)
     if(nbend==1) then
        j1 = i - 1
        j2 = i
        j3 = i + 1
        rxyz1 = (/h_X(j1),h_Y(j1),h_Z(j1)/)
        rxyz2 = (/h_X(j2),h_Y(j2),h_Z(j2)/)
        rxyz3 = (/h_X(j3),h_Y(j3),h_Z(j3)/)
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
        beta(i) = beta_tmp
        !write(*,*) 'beta_tmp',beta_tmp,cosbet,dotpr
     endif
  enddo
  
  return
end subroutine updateLH

