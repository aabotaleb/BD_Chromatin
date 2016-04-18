
      subroutine readcore(corefile,Nq,Nq3,core_pos,core_q,myid)

      implicit NONE

      integer Nq, Nq3,myid
      double precision core_pos(Nq3), core_q(Nq)
      character  corefile*50

      integer i, i1,i2,i3
c    salt concentration dependent core_data
      open (unit = 1, file = corefile, 
     + form = 'formatted', access = 'sequential', status = 'old')
     
       write(*,'(A,i2)'),' COREFILE BEING READ BY........',myid

      do 50 i = 1,Nq
        i1 = 3*(i-1)+1
        i2 = i1+1
        i3 = i2+1
        read (unit = 1, fmt = 40, end = 100)
     +        core_pos(i1), core_pos(i2), core_pos(i3)
   40   format( e14.7, 1x, e14.7, 1x, e14.7, 1x )
c ..... convert Angstrom to nm: ...........
        core_pos(i1) = core_pos(i1) / 10.0d0
        core_pos(i2) = core_pos(i2) / 10.0d0
        core_pos(i3) = core_pos(i3) / 10.0d0
c ..........................................
   50 continue

      do 90 i = 1,Nq
        read (unit = 1, fmt = 80, end = 100) core_q(i)
   80   format( e14.7)
   90 continue

      close ( unit = 1 )
 
      return

100   print*, 'Premature end of file in subroutine: readcore'

      end







