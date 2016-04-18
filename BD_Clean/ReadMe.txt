
Analytical and Numerical Calculations of the Mechanical forces 
for the chromatin only are done 
[All mechanical forces : Stretching - Bending and Twistting]
Plus testing the correctness of the results along with testgh algorithm[Based on Taylor expansion]

1- type : make 
to run the makefile and so compiling the fortran files 
and generating the executable : "chrom_vNRL.x"

2- type ./run.sh to run the script "run.sh" 

3-the output files after run are : 
	1-output = "log file"
	2-Analytically calculated forces = "Analytic Forces[MATLAB Version].txt"
	3-Numerically calculated forces  = "Numerical Forces.txt"
	4-testgh algorithm output results to check the correctness of the numerical differentation:
		i-Stretching mechanical forces  : "TestRes_Stretching_Numerical.txt"
		ii-Bending mechanical forces    : "TestRes_Bending_Numerical.txt"
		iii-Twistting mechanical forces : "TestRes_Twisting_Numerical.txt"
	If RATIO tends to 4 or 8 as eps is decreased,G  is correct or G&H are correct
	where G is the gradient and H is the hessian and as we are conerned with the potential Gradient
	so when RATIO tends to 4 it confirms that the gradient is correct.
	5-Coordiantes and euler angles of the core and linker beads : "Coordinates data xx . txt " where xx is the bead number.
	6-Energies of the core and linker beads : "Energies data  xx .txt " where xx is the bead number.
		