cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine startconf(n_c,nlb,n,n3,type,r,ro,a,b,c,lo,np,myid,ierr)

C      use mpi
      implicit NONE

c     beads of core and linker DNA on each repeating unit.
      integer n_c,n,n3, type(n),nlb(n_c)
      double precision r(n3),ro, a(n3),b(n3),c(n3),lo
      double precision norma,normb,normc
      integer i, j, j1, j2, j3, index,m1,m2,m3
      integer j1_c, j2_c, j3_c

      double precision ac(3),bc(3),cc(3),ccc(3)
      double precision ac_old(3),bc_old(3),cc_old(3)
      double precision norm,proj
      integer i1,i2,i3

      integer np,myid,ierr
	  
C	  print*,'entering start_conf on proc',myid
   

c      1.dat contains position and a, b, c of nucleosome cores and linkers
      open (unit = 22, file = '10.dat', form = 'formatted',
     +      access = 'sequential', status = 'old')
      write(*,'(A,i2)'),' COORDINATE FILE BEING READ BY.',myid
c     read linker bead coordinates (***6-bead linkers only)
      index=0    
 


      index=0
      do 10 i = 1, n_c
         index=index+1
         type(index)=1
         if(i.eq.n_c)type(index)=2
       do 15 j=1,nlb(i)
            j1 = 3*index + 1
            j2 = 3*index + 2
            j3 = 3*index + 3
            read(22,*) r(j1), r(j2), r(j3)
            index=index+1
            type(index)=0
 15      continue
 10   continue



      index=0
      do 20 i = 1, n_c
         j1 = 3*index + 1
         j2 = 3*index + 2
         j3 = 3*index + 3
        read(22,*) r(j1), r(j2), r(j3)
   	  read(22,*) a(j1), a(j2), a(j3)
	  read(22,*) b(j1), b(j2), b(j3)
	  read(22,*) c(j1), c(j2), c(j3)
          norma=dsqrt(a(j1)**2+a(j2)**2+a(j3)**2)
          normb=dsqrt(b(j1)**2+b(j2)**2+b(j3)**2)
          normc=dsqrt(c(j1)**2+c(j2)**2+c(j3)**2)
          a(j1)=a(j1)/norma         
          a(j2)=a(j2)/norma         
          a(j3)=a(j3)/norma         

          b(j1)=b(j1)/normb         
          b(j2)=b(j2)/normb         
          b(j3)=b(j3)/normb        
 
          c(j1)=c(j1)/normc         
          c(j2)=c(j2)/normc         
          c(j3)=c(j3)/normc        
 
         index=index+nlb(i)+1

 20   continue


      index=1
      do i = 1, n
         if(type(i).eq.0)then
             j1 = 3*(i-1) + 1
             j2 = 3*(i-1) + 2
             j3 = 3*(i-1) + 3
	     a(j1) = a(3*(index-1) + 1)
	     a(j2) = a(3*(index-1) + 2)
	     a(j3) = a(3*(index-1) + 3)
	     b(j1) = b(3*(index-1) + 1)
	     b(j2) = b(3*(index-1) + 2)
	     b(j3) = b(3*(index-1) + 3)
	     c(j1) = c(3*(index-1) + 1)
	     c(j2) = c(3*(index-1) + 2)
	     c(j3) = c(3*(index-1) + 3)
	 else
	     index = i
	 endif   
      enddo




ccccccccccc
ccc Modification: DURBA and TONI
c     We make sure that a,b,c are orthonormal 
      do i = 1,n
c      do i = 1,1
         i1 = 3*(i-1) + 1
         i2 = i1 + 1
         i3 = i2 + 1

         ac = (/a(i1),a(i2),a(i3)/)
         bc = (/b(i1),b(i2),b(i3)/)
         cc = (/c(i1),c(i2),c(i3)/)

         ac_old = ac
         bc_old = bc
         cc_old = cc

c     Renormalize 'a'
         norm = sqrt(dot_product(ac,ac))
         ac = ac/norm
c     Correct 'b' and renormalize
         proj = dot_product(ac,bc)
         bc = bc - proj*ac
         norm = sqrt(dot_product(bc,bc))
         bc = bc/norm
c     Cross product for 'c' and renormalize
         cc(1) = ac(2)*bc(3) - ac(3)*bc(2)
         cc(2) = ac(3)*bc(1) - ac(1)*bc(3)
         cc(3) = ac(1)*bc(2) - ac(2)*bc(1)
         norm = sqrt(dot_product(cc,cc))
         cc = cc/norm	            

         a(i1) = ac(1)
         a(i2) = ac(2)
         a(i3) = ac(3)

         b(i1) = bc(1)
         b(i2) = bc(2)
         b(i3) = bc(3)

         c(i1) = cc(1)
         c(i2) = cc(2)
         c(i3) = cc(3)


      end do
cccccccccc
   
C      print*, 'AFTER STARTCONF'
C      write(*,*) 'proc,r=',myid,r(1)
C	  write(*,*) 'a after start_conf=',myid,a(1)
C	  write(*,*) 'b=',myid,b(1)
C	  write(*,*) 'c=',myid,c(1)

c
      return
100   print*, 'Premature end of file in subroutine: build.f'
      end
