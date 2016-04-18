### MAKEFILE for chromatin code on MU-FE
### MPI libraries given below, for performance tips see README
### Mustang@LANL on Turquoise Network 




###MOD = /usr/projects/hpcsoft/toss2/common/intel-clusterstudio/2015.0.0.022/impi_5.0.1/bin64/mpif77 -c modglob.f
###CC= /usr/projects/hpcsoft/toss2/common/intel-clusterstudio/2015.0.0.022/impi_5.0.1/bin64/mpif77 -I/usr/projects/hpcsoft/toss2/common/intel-clusterstudio/2015.0.0.022/impi_5.0.1/include/ -L//usr/projects/hpcsoft/toss2/common/intel-clusterstudio/2015.0.0.022/impi_5.0.1/lib/ 
RM1 = rm -f modglob.mod
RM2 = rm -f chrom_vNRL.x
MOD = mpif90 -c modglob.f
CC= mpif90 -g -O3 -pg 


SRC =   main.f potential.f startconf.f readcore.f readtail.f \
        update_mod.f mersenne.f \
        modglob.o \
        LHboundsub.f90  testgh.f\
	gfat.f init.f 

OBJS =  main.o potential.o startconf.o  readcore.o readtail.o \
        update_mod.o  mersenne.o\
        modglob.o \
        LHboundsub.o testgh.o  \
	gfat.o init.o 

compile:
	$(RM1)
	$(RM2)
	$(MOD)
	$(CC) $(SRC)  -o chrom_vNRL.x
clean:
	rm -f $(OBJS) 

veryclean:
	 rm *.o bd
