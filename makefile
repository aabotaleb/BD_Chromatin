### MAKEFILE for chromatin code on MU-FE
### MPI libraries given below, for performance tips see README
### Mustang@LANL on Turquoise Network 




###MOD = /usr/projects/hpcsoft/toss2/common/intel-clusterstudio/2015.0.0.022/impi_5.0.1/bin64/mpif77 -c modglob.f
###CC= /usr/projects/hpcsoft/toss2/common/intel-clusterstudio/2015.0.0.022/impi_5.0.1/bin64/mpif77 -I/usr/projects/hpcsoft/toss2/common/intel-clusterstudio/2015.0.0.022/impi_5.0.1/include/ -L//usr/projects/hpcsoft/toss2/common/intel-clusterstudio/2015.0.0.022/impi_5.0.1/lib/ 
RM1 = rm -f modglob.mod
RM2 = rm -f chrom_vNRL.x
MOD = mpif90 -c modglob.f
CC= mpif90 -g -O3 -pg -heap-arrays


SRC =   main.f potential.f startconf.f readcore.f readtail.f \
        update_mod.f global.f  mersenne.f \
        grids.f \
        modglob.o LHbead.f \
        updateLH.f90 LHboundsub.f90  testgh.f calcfat.f \
		diffcent.f  geteuler.f gfat.f init.f Isimpson.f

OBJS =  main.o potential.o startconf.o  readcore.o readtail.o \
        update_mod.o global.o mersenne.o\
        grids.o \
        modglob.o LHbead.o \
        updateLH.o LHboundsub.o testgh.o calcfat.o \
		diffcent.o geteuler.o gfat.o init.o Isimpson.o


compile:
	$(RM1)
	$(RM2)
	$(MOD)
	$(CC) $(SRC)  -o chrom_vNRL.x
clean:
	rm -f $(OBJS) 

veryclean:
	 rm *.o bd
