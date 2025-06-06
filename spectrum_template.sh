#!/bin/sh -l

#SBATCH -A mth240048
#SBATCH -p wholenode
#SBATCH -o goldfish.o%j     # Name of stdout output file
#SBATCH -e goldfish.e%j     # Name of stderr error file
#SBATCH --mail-user=liamoconnor2025@u.northwestern.edu
#SBATCH --mail-type=all   # Send email to above address at begin and end of job



source /home/x-loconnor/.bashrc
conda activate dedalus3
export MPI_UNBUFFERED_STDIO=true

source ~/png2mp4.sh

FILE="$(readlink -f "$0")"
DIR="$(dirname "$(readlink -f "$0")")/"
CONFIG="options.cfg"
# SOLVER="rpcf-mhd.py"


# echo "RESOURCES ALLOCATED: $RLARG" > message.txt
# echo "RESOURCES ALLOCATED: $RLARG" 
echo "OG SUFIX SUPPLIED: $suffix"
echo "currently in (skip. her's ls)"
echo $DIR
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo $SCRIPT_DIR

mkdir traces
ls
cat options.cfg > /home/x-loconnor/temp.cfg
source /home/x-loconnor/temp.cfg
rm /home/x-loconnor/temp.cfg
echo "CONFIG SUFIX SUPPLIED: $suffix" 

echo "ky, growth_rate, std_err" > root_output0.txt
for num in $(seq 0.005 0.005 1.0); do
    # echo $num
    mpirun -np $MPIPROC python3 $SOLVER $num
    # exit 1
done