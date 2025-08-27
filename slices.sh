#!/bin/bash
# module load ffmpeg
source ~/png2mp4.sh

function make_video() {
    # png2mp4_big $source $videoname $fps
    png2mp4 $source $videoname $fps
    echo "~/tvd-fast/$videoname"
}

suffix=$1
# procs=128
procs=$2
# export LD_PRELOAD=/usr/lib64/libslurm.so
# mpirun -np 1 python3 plot_energies.py $suffix
# fps=360
echo $procs

fps=120
mpirun -np $procs python3 plot_clean.py $suffix
mpirun -np $procs python3 plot_slicepoint.py $suffix


source="$suffix/cleany/"
videoname="$suffix/cleany.mp4"
make_video

source="$suffix/mid_zx/"
videoname="$suffix/mid_zx.mp4"
make_video

source="$suffix/cleanx/"
videoname="$suffix/cleanx.mp4"
make_video

source="$suffix/cleanz/"
videoname="$suffix/cleanz.mp4"
make_video
exit 1
# exit 1

source="$suffix/mid_yz/"
videoname="$suffix/mid_yz.mp4"
make_video

source="$suffix/mid_xy/"
videoname="$suffix/mid_xy.mp4"
make_video


source="$suffix/avg_xy/"
videoname="$suffix/avg_xy.mp4"
make_video

source="$suffix/avg_yz/"
videoname="$suffix/avg_yz.mp4"
make_video


source="$suffix/avg_zx/"
videoname="$suffix/avg_zx.mp4"
make_video

source="$suffix/profiles_avg/"
videoname="$suffix/profles_avg.mp4"
make_video

exit 1


# png2mp4 frames_yz/ midplane_yz.mp4 $fps
# png2mp4 frames_zx/ midplane_zx.mp4 $fps

