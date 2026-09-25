#!/bin/sh -l
#SBATCH -A mth240048
#SBATCH -p wholenode
#SBATCH -o goldfish.o%j     # Name of stdout output file
#SBATCH -e goldfish.e%j     # Name of stderr error file
#SBATCH --mail-user=liamoconnor2025@u.northwestern.edu
#SBATCH --mail-type=all   # Send email to above address at begin and end of job
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=0:30:00
#SBATCH --job-name slices_combined

# Submit with the run-directory suffixes you'd otherwise pass to slices.sh, e.g.:
#   sbatch submit_slices_combined.sh Ro3p5_Lz8pi_Ny16_Rm175 Ro3p5_Lz8pi_Ny16_Rm175_RSTRT1
# Everything after the script name is forwarded to slices.sh as "$@".

source /home/x-loconnor/.bashrc
conda activate dedalus3
export MPI_UNBUFFERED_STDIO=true

cd ~/tvd-fast/
bash slices.sh "$@"
