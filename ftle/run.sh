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
# DIR="$(dirname "$(readlink -f "$0")")/"

PROCS=16
cd ~/tvd-fast/ftle

# mpirun -n $PROCS python3 shear_paths_new.py --restartN=72 --name=NW_72
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=72 --name=NE_72
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=72 --name=SW_72
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=72 --name=SE_72

# mpirun -n $PROCS python3 shear_paths_new.py --restartN=40 --name=NW_40
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=40 --name=NE_40
# mpirun -n $PROCS python3 shear_paths_new.py --restartN=40 --name=SW_40
mpirun -n $PROCS python3 shear_paths_new.py --restartN=40 --name=SE_40