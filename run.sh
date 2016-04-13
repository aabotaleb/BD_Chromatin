##############################
#!/bin/bash
#SBATCH --job-name=test_chromatin
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --ntasks=24
#SBATCH --time=20:00

source /etc/profile.d/modules.sh
module load intel/12.0.0
module load openmpi/intel/1.10.2

cd /home/mamin/BD/BD_Only
srun --ntasks-per-node=12 mpirun -np 2 ./chrom_vNRL.x > output

#################################
