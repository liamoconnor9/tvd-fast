#!/bin/bash
# module load ffmpeg
source ~/png2mp4.sh
source /home/x-loconnor/.bashrc
conda activate dedalus3
# exit 1
function make_video() {
    # png2mp4_big $source $videoname $fps
    png2mp4 $source $videoname $fps
    echo "~/tvd-fast/$videoname"
}

# Usage:
#   bash slices.sh <suffix> [procs] [--stride=N] [--start-time=T] [--end-time=T] [--dark] [--kezmodes] [--energies]
#   bash slices.sh <suffix1> <suffix2> [...] [procs] [--stride=N] [--start-time=T] [--end-time=T] [--dark] [--kezmodes] [--energies]
#   bash slices.sh <suffix> --chain [procs] [--stride=N] [--start-time=T] [--end-time=T] [--dark] [--kezmodes] [--energies]
#                                                          (auto-discovers
#                                                           <suffix>_RSTRT1, _RSTRT2, ...)
#
# --dark renders with plot_clean.py's dark-background style (see plot_clean.py
# --dark); omitted, frames use the default light/clean style.
#
# --start-time=T / --end-time=T only plot/use frames whose simulation time
# (the same "t = ..." value burned into each frame, i.e. continuous across
# restarts in combined mode) falls in [T_start, T_end]; either bound can be
# given alone. Like --stride, frames outside the range are simply never
# rendered on a fresh run, but if frames were already rendered in a previous
# run, selection for merging/compositing is still done from what's already
# there, so no re-render is needed just to preview a narrower time window.
#
# In combined mode (2+ dirs, or --chain), cleany/cleanx frames are
# regenerated per directory as usual, then stitched into one
# chronologically-ordered movie per plane (trimming each run but the last to
# its own final checkpoint time, since every restart in this codebase resets
# its write-number/time counters to zero -- see merge_slice_frames.py) and
# saved as {plane}_combined.mp4 in the first directory only (the original
# simulation). The "t = ..." burned into each frame is also shifted per
# directory (merge_slice_frames.py --offsets) so the time shown in the
# combined movie is continuous across restarts instead of resetting to 0 at
# each seam.
#
# --stride=N only plots/uses every Nth frame (write_number 1, 1+N, 1+2N, ...),
# for a faster/cheaper preview video. Frames not selected by the stride are
# simply never rendered (see plot_clean.py --stride), so this saves real
# render time on a fresh run -- but it also works if every frame was already
# rendered in a previous (stride=1 or different-stride) run: frame selection
# for the movie is done by filename (write_number), not by what got rendered
# this time, so existing frames are just subsampled, nothing is deleted or
# needs to be re-rendered.
#
# --kezmodes additionally renders {plane}_kezmodes.mp4 (or
# {plane}_combined_kezmodes.mp4 in combined mode): the run's kezmodes plot
# (ke_mode1..10 vs time, spanning the whole chain) shown beside the normal
# slice plot, with a vertical line sweeping across as the video plays to mark
# the current frame's time. Both panels are rendered into one matplotlib
# figure per frame (see plot_kezmodes_combined.py) rather than pasting two
# separately-rendered images together, so the kezmodes panel's frame lines up
# exactly with the slice panel's frame and the marker line is bounded to it.
# Saved in the first directory only, same as the plain combined video.
#
# --energies behaves exactly like --kezmodes (same per-frame single-figure
# rendering, same "saved in the first directory only" placement, and it also
# suppresses the plain video, same as --kezmodes) but the side panel instead
# shows the run's energies plot (be_x/y/z, ke_x/y/z vs time, log-scaled --
# same six tasks/colors as energies__combined.png). Renders
# {plane}_energies.mp4 (or {plane}_combined_energies.mp4 in combined mode).
#
# --kezmodes and --energies together render ONE video (not two): the two
# panels are stacked vertically in the same right-hand column (energies on
# top, kezmodes below, sharing one x-axis) -- see plot_kezmodes_combined.py.
# Saved as {plane}_kezmodes_energies.mp4 (or
# {plane}_combined_kezmodes_energies.mp4 in combined mode).

DIRS=()
PROCS=""
CHAIN=false
STRIDE=1
DARK=false
KEZMODES=false
ENERGIES=false
START_TIME=""
END_TIME=""
for arg in "$@"; do
    if [ "$arg" = "--chain" ] || [ "$arg" = "-c" ]; then
        CHAIN=true
    elif [ "$arg" = "--dark" ] || [ "$arg" = "-d" ]; then
        DARK=true
    elif [ "$arg" = "--kezmodes" ] || [ "$arg" = "-k" ]; then
        KEZMODES=true
    elif [ "$arg" = "--energies" ] || [ "$arg" = "-e" ]; then
        ENERGIES=true
    elif [[ "$arg" == --stride=* ]]; then
        STRIDE="${arg#--stride=}"
    elif [[ "$arg" == --start-time=* ]]; then
        START_TIME="${arg#--start-time=}"
    elif [[ "$arg" == --end-time=* ]]; then
        END_TIME="${arg#--end-time=}"
    elif [[ "$arg" =~ ^[0-9]+$ ]]; then
        PROCS="$arg"
    else
        DIRS+=("$arg")
    fi
done
# Default procs to whatever's actually available right now: inside a Slurm
# allocation, use the cores it granted; otherwise (e.g. running loose on a
# login/interactive node) fall back to the node's logical core count.
if [ -z "$PROCS" ]; then
    if [ -n "$SLURM_NTASKS" ]; then
        PROCS="$SLURM_NTASKS"
    elif [ -n "$SLURM_CPUS_ON_NODE" ]; then
        PROCS="$SLURM_CPUS_ON_NODE"
    else
        PROCS=$(nproc)
    fi
fi
fps=60

DARK_FLAG=""
if $DARK; then
    DARK_FLAG="--dark"
fi

START_TIME_FLAG=""
if [ -n "$START_TIME" ]; then
    START_TIME_FLAG="--start-time=$START_TIME"
fi

END_TIME_FLAG=""
if [ -n "$END_TIME" ]; then
    END_TIME_FLAG="--end-time=$END_TIME"
fi

STRIDE_SUFFIX=""
if [ "$STRIDE" -gt 1 ]; then
    STRIDE_SUFFIX="_stride${STRIDE}"
fi

TIME_SUFFIX=""
if [ -n "$START_TIME" ] || [ -n "$END_TIME" ]; then
    TIME_SUFFIX="_t${START_TIME:-0}-${END_TIME:-end}"
fi

BOTH=false
if $KEZMODES && $ENERGIES; then
    BOTH=true
fi

if $CHAIN; then
    if [ "${#DIRS[@]}" -ne 1 ]; then
        echo "--chain expects exactly one suffix to extend"
        exit 1
    fi
    base="${DIRS[0]%/}"
    n=1
    while [ -d "${base}_RSTRT${n}" ]; do
        DIRS+=("${base}_RSTRT${n}")
        n=$((n+1))
    done
    if [ "${#DIRS[@]}" -eq 1 ]; then
        echo "no ${base}_RSTRT* restarts found; nothing to chain"
        exit 1
    fi
    echo "auto-discovered restart chain: ${DIRS[@]}"
fi

if [ "${#DIRS[@]}" -eq 1 ]; then
    suffix="${DIRS[0]}"
    procs=$PROCS
    echo $procs

    mpirun -np $procs python3 plot_clean.py $suffix $DARK_FLAG --stride=$STRIDE $START_TIME_FLAG $END_TIME_FLAG

    if [ "$STRIDE" -gt 1 ] || [ -n "$START_TIME_FLAG" ] || [ -n "$END_TIME_FLAG" ] || $KEZMODES || $ENERGIES; then
        for plane in cleany cleanx; do
            if ! $KEZMODES && ! $ENERGIES; then
                staging=$(mktemp -d)
                python3 merge_slice_frames.py $plane $staging --stride=$STRIDE $START_TIME_FLAG $END_TIME_FLAG "$suffix"
                if [ $? -eq 0 ]; then
                    source="$staging/"
                    videoname="$suffix/${plane}${STRIDE_SUFFIX}${TIME_SUFFIX}.mp4"
                    make_video
                else
                    echo "skipping $plane: no frames to combine"
                fi
                rm -rf $staging
            fi

            if $BOTH; then
                staging4=$(mktemp -d)
                python3 plot_kezmodes_combined.py $plane $staging4 --stride=$STRIDE $DARK_FLAG --kezmodes --energies $START_TIME_FLAG $END_TIME_FLAG "$suffix"
                if [ $? -eq 0 ]; then
                    source="$staging4/"
                    videoname="$suffix/${plane}_kezmodes_energies${STRIDE_SUFFIX}${TIME_SUFFIX}.mp4"
                    make_video
                else
                    echo "skipping $plane kezmodes+energies composite: no frames"
                fi
                rm -rf $staging4
            else
                if $KEZMODES; then
                    staging2=$(mktemp -d)
                    python3 plot_kezmodes_combined.py $plane $staging2 --stride=$STRIDE $DARK_FLAG $START_TIME_FLAG $END_TIME_FLAG "$suffix"
                    if [ $? -eq 0 ]; then
                        source="$staging2/"
                        videoname="$suffix/${plane}_kezmodes${STRIDE_SUFFIX}${TIME_SUFFIX}.mp4"
                        make_video
                    else
                        echo "skipping $plane kezmodes composite: no frames"
                    fi
                    rm -rf $staging2
                fi

                if $ENERGIES; then
                    staging3=$(mktemp -d)
                    python3 plot_kezmodes_combined.py $plane $staging3 --stride=$STRIDE $DARK_FLAG --energies $START_TIME_FLAG $END_TIME_FLAG "$suffix"
                    if [ $? -eq 0 ]; then
                        source="$staging3/"
                        videoname="$suffix/${plane}_energies${STRIDE_SUFFIX}${TIME_SUFFIX}.mp4"
                        make_video
                    else
                        echo "skipping $plane energies composite: no frames"
                    fi
                    rm -rf $staging3
                fi
            fi
        done
        exit 0
    fi

    source="$suffix/cleany/"
    videoname="$suffix/cleany.mp4"
    make_video

    source="$suffix/cleanx/"
    videoname="$suffix/cleanx.mp4"
    make_video

    source="$suffix/cleanz/"
    videoname="$suffix/cleanz.mp4"
    make_video
    exit 0
fi

if [ "${#DIRS[@]}" -eq 0 ]; then
    echo "usage: bash slices.sh <suffix1> [<suffix2> ...] [procs] [--stride=N] [--start-time=T] [--end-time=T] [--dark] [--kezmodes] [--energies]"
    echo "       bash slices.sh <suffix> --chain [procs] [--stride=N] [--start-time=T] [--end-time=T] [--dark] [--kezmodes] [--energies]"
    exit 1
fi

echo "combined mode, procs=$PROCS, dirs=${DIRS[@]}, stride=$STRIDE"

mapfile -t OFFSETS < <(python3 merge_slice_frames.py --offsets "${DIRS[@]}")

for i in "${!DIRS[@]}"; do
    d="${DIRS[$i]}"
    off="${OFFSETS[$i]}"
    mpirun -np $PROCS python3 plot_clean.py $d $DARK_FLAG --time-offset=$off --stride=$STRIDE $START_TIME_FLAG $END_TIME_FLAG
done

for plane in cleany cleanx; do
    if ! $KEZMODES && ! $ENERGIES; then
        staging=$(mktemp -d)
        python3 merge_slice_frames.py $plane $staging --stride=$STRIDE $START_TIME_FLAG $END_TIME_FLAG "${DIRS[@]}"
        if [ $? -eq 0 ]; then
            source="$staging/"
            videoname="${DIRS[0]}/${plane}_combined${STRIDE_SUFFIX}${TIME_SUFFIX}.mp4"
            png2mp4 $source $videoname $fps
            echo "$videoname"
        else
            echo "skipping $plane: no frames to combine"
        fi
        rm -rf $staging
    fi

    if $BOTH; then
        staging4=$(mktemp -d)
        python3 plot_kezmodes_combined.py $plane $staging4 --stride=$STRIDE $DARK_FLAG --kezmodes --energies $START_TIME_FLAG $END_TIME_FLAG "${DIRS[@]}"
        if [ $? -eq 0 ]; then
            source="$staging4/"
            videoname="${DIRS[0]}/${plane}_combined_kezmodes_energies${STRIDE_SUFFIX}${TIME_SUFFIX}.mp4"
            png2mp4 $source $videoname $fps
            echo "$videoname"
        else
            echo "skipping $plane kezmodes+energies composite: no frames"
        fi
        rm -rf $staging4
    else
        if $KEZMODES; then
            staging2=$(mktemp -d)
            python3 plot_kezmodes_combined.py $plane $staging2 --stride=$STRIDE $DARK_FLAG $START_TIME_FLAG $END_TIME_FLAG "${DIRS[@]}"
            if [ $? -eq 0 ]; then
                source="$staging2/"
                videoname="${DIRS[0]}/${plane}_combined_kezmodes${STRIDE_SUFFIX}${TIME_SUFFIX}.mp4"
                png2mp4 $source $videoname $fps
                echo "$videoname"
            else
                echo "skipping $plane kezmodes composite: no frames"
            fi
            rm -rf $staging2
        fi

        if $ENERGIES; then
            staging3=$(mktemp -d)
            python3 plot_kezmodes_combined.py $plane $staging3 --stride=$STRIDE $DARK_FLAG --energies $START_TIME_FLAG $END_TIME_FLAG "${DIRS[@]}"
            if [ $? -eq 0 ]; then
                source="$staging3/"
                videoname="${DIRS[0]}/${plane}_combined_energies${STRIDE_SUFFIX}${TIME_SUFFIX}.mp4"
                png2mp4 $source $videoname $fps
                echo "$videoname"
            else
                echo "skipping $plane energies composite: no frames"
            fi
            rm -rf $staging3
        fi
    fi
done
