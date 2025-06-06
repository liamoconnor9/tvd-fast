#!/bin/sh -l

#SBATCH --job-name rm1500ky1 
#SBATCH -A mth240048
#SBATCH -p wholenode
#SBATCH -o goldfish.o%j     # Name of stdout output file
#SBATCH -e goldfish.e%j     # Name of stderr error file
#SBATCH --mail-user=liamoconnor2025@u.northwestern.edu
#SBATCH --mail-type=all   # Send email to above address at begin and end of job

# shopt -s expand_aliases
# alias mpiexec_mpt="mpirun"
# alias ffmpeg3="ffmpeg"

# export PATH=$HOME/scripts:$PATH
# deactivate
# unset PYTHONPATH
# export PYTHONNOUSERSITE=1
# support lots of text output to stdio for analysis

# wholenode
# highmem

source /home/x-loconnor/.bashrc
conda activate dedalus3
export MPI_UNBUFFERED_STDIO=true

source ~/png2mp4.sh

FILE="$(readlink -f "$0")"
DIR="$(dirname "$(readlink -f "$0")")/"
CONFIG="options.cfg"
SOLVER="rpcf-mhd.py"


# echo "RESOURCES ALLOCATED: $RLARG" > message.txt
# echo "RESOURCES ALLOCATED: $RLARG" 
echo "OG SUFIX SUPPLIED: $suffix"
echo "currently in (skip. her's ls)"
echo $DIR
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo $SCRIPT_DIR

ls
cat options.cfg > /home/x-loconnor/temp.cfg
source /home/x-loconnor/temp.cfg
rm /home/x-loconnor/temp.cfg
echo "CONFIG SUFIX SUPPLIED: $suffix" 

# for seed in $(seq $sim_index $Nsims)
# do 
# echo "##############################################################"
# # echo "              Running Simulation With Seed $seed"
# echo "##############################################################"
mpirun -np $MPIPROC python3 $SOLVER $CONFIG
# cd ..
# exit 1
# bash slices.sh $suffix $MPIPROC
# # mpiexec_mpt -np $MPIPROC python3 plot_slicepoint.py $suffix
# # mpiexec_mpt -np 1 python3 plot_energies.py $suffix

# done


# # sed -i 's/seed=42/seed=2/' ~/mhd/mri_options.cfg

# # mpiexec_mpt -np $MPIPROC python3 mri.py $CONFIG