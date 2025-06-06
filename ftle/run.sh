#!/bin/sh -l
#SBATCH -A mth240048
#SBATCH -p wholenode
#SBATCH -o goldfish.o%j     # Name of stdout output file
#SBATCH -e goldfish.e%j     # Name of stderr error file
#SBATCH --mail-user=liamoconnor2025@u.northwestern.edu
#SBATCH --mail-type=all   # Send email to above address at begin and end of job
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=96:00:00
#SBATCH --job-name quadrants_test

source /home/x-loconnor/.bashrc
conda activate dedalus2
export MPI_UNBUFFERED_STDIO=true

source ~/png2mp4.sh

FILE="$(readlink -f "$0")"
DIR="$(dirname "$(readlink -f "$0")")/"

PROCS=16
cd $DIR

# mpirun -n $PROCS python3 shear_paths_new.py --restartN=33 --name=NW_33
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=33 --name=NE_33
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=33 --name=SW_33
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=33 --name=SE_33

# mpirun -n $PROCS python3 shear_paths_new.py --restartN=1 --name=NW_1
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=1 --name=NE_1
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=1 --name=SW_1
mpirun -n $PROCS python3 shear_paths_new.py --restartN=1 --name=SE_1