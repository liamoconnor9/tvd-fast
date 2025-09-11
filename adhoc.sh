#!/bin/sh -l
#SBATCH -A mth240048
#SBATCH -p wholenode
#SBATCH -o goldfish.o%j     # Name of stdout output file
#SBATCH -e goldfish.e%j     # Name of stderr error file
#SBATCH --mail-user=liamoconnor2025@u.northwestern.edu
#SBATCH --mail-type=all   # Send email to above address at begin and end of job
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=4:00:00
#SBATCH --job-name plot

source /home/x-loconnor/.bashrc
conda activate dedalus3
export MPI_UNBUFFERED_STDIO=true

source ~/png2mp4.sh

FILE="$(readlink -f "$0")"
DIR="$(dirname "$(readlink -f "$0")")/"

PROCS=128
# cd $DIR
cd ~/tvd-fast/
bash slices.sh nl_Lz8pi_Rm3232 $PROCS