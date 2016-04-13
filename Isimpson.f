Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
      subroutine   Isimpson(n,t_X,Y,res)   

      implicit none  !CSUN 6/2/04
      integer i,n
      double precision Y(n)
      double precision t_X(n)
      double precision res
      res=0.d0
     
      ! trapeziodal unequal spaced 
       do i=1,n-1
           res=res+((t_X(i+1)-t_X(i))*0.50d0*(Y(i)+Y(i+1)))
       enddo
      end
