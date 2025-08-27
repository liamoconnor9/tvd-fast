#!/bin/sh -l

#SBATCH --job-name test-standard 
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

source ~/.bashrc
conda activate dedalus3
export MPI_UNBUFFERED_STDIO=true

source ~/png2mp4.sh

FILE="$(readlink -f "$0")"
DIR="$(dirname "$(readlink -f "$0")")/"
CONFIG="options.cfg"
echo "OG SUFIX SUPPLIED: $suffix"
echo "currently in (skip. her's ls)"
echo $DIR
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo $SCRIPT_DIR
cat options.cfg > ~/lksdjfsl.cfg
source ~/lksdjfsl.cfg
rm ~/lksdjfsl.cfg
echo "CONFIG SUFIX SUPPLIED: $suffix" 
mpirun -np $MPIPROC python3 $SOLVER $CONFIG
