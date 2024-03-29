	subroutine grids(Ngrids,gridsize,boundx,boundy,boundz,occmatrix,myid)
	
	implicit NONE
	
	integer			i,j,k,Ngrids, occmatrix(14,14,10)
	integer			locx,locy,locz,myid
	double precision	gridsize,boundx,boundy,boundz
	

	do i = 1, 14
	   do j = 1, 14
	      do k = 1, 10
	         occmatrix(i,j,k) = 0
	      enddo
	   enddo
	enddo
	
	open (13, file = 'grid_data.9deg')	
		 write(*,'(A,i2)'),' READING GRID FILE ON PROC.....',myid
        read(13,*) Ngrids,gridsize
	read(13,*) boundx,boundy,boundz		
	do i = 1, Ngrids
	    read(13,*) locx, locy, locz
	    occmatrix(locx,locy,locz) = 1 
	enddo	
	close(13)
	
	return
	end
