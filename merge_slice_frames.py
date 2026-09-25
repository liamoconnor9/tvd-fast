import sys
import os
import re
import glob
import h5py

# Builds an ordered, de-duplicated sequence of slice-plot frames (as produced
# by plot_clean.py) across a chain of restarted runs, and symlinks them into a
# staging directory with fresh sequential filenames so that a plain
# `cat staging_dir/* | ffmpeg ...` (see png2mp4 in ~/png2mp4.sh) encodes them
# in the correct chronological order.
#
# Every run in this codebase resets both solver.sim_time and the analysis
# write_number counter to zero on restart (see rpcf-mhd.py / the
# plot_energies_combined.py writeup), so frame numbers collide across run
# directories and don't reflect true simulation time. To combine:
#   - every run except the last is trimmed to the sim_time of its own final
#     checkpoint (the point the next run in the chain actually restarted
#     from -- frames logged after that were never resumed from)
#   - frames are then taken in run order, and within each run in increasing
#     write_number order (write_number increases monotonically with
#     sim_time within a single run)
#
# Usage:
#   python3 merge_slice_frames.py <plane_tag> <staging_dir> <dir1> <dir2> [<dir3> ...]
# e.g.
#   python3 merge_slice_frames.py cleany /tmp/stage123 Ro3p5_Lz8pi_Ny16_Rm175 Ro3p5_Lz8pi_Ny16_Rm175_RSTRT1

def get_checkpoint_final_time(dir):
    """sim_time of the last (highest-index) checkpoint saved in `dir`."""
    files = glob.glob("{}checkpoint/checkpoint_s*.h5".format(dir))
    if not files:
        return None
    def idx(f):
        return int(os.path.basename(f).split('checkpoint_s')[-1].split('.h5')[0])
    files.sort(key=idx)
    with h5py.File(files[-1], "r") as f:
        return float(f['scales']['sim_time'][()][0])

def natural_sorted_h5(pattern):
    files = glob.glob(pattern)
    def idx(f):
        m = re.search(r'_s(\d+)\.h5$', os.path.basename(f))
        return int(m.group(1))
    files.sort(key=idx)
    return files

def frame_list(dir, cutoff, stride=1, offset=0.0, start_time=None, end_time=None):
    """[(write_number, sim_time), ...] for slicepoints writes in `dir`, in
    chronological order, trimmed to sim_time <= cutoff when cutoff is given,
    and subsampled to every `stride`-th write (write_number 1, 1+stride, ...).
    This is the same predicate plot_clean.py's --stride uses to decide which
    frames to render, so a directory that was only ever rendered with a given
    stride has exactly the frames this selects -- and a directory rendered
    densely (every write) works too, since this just subsamples what's there.

    start_time/end_time bound the *combined/global* time (st + offset, same
    as the "t = ..." burned into each frame by plot_clean.py's
    --time-offset), so they mean the same thing here as in plot_clean.py
    regardless of which run in a restart chain a frame actually came from."""
    entries = []
    for file in natural_sorted_h5("{}slicepoints/*.h5".format(dir)):
        with h5py.File(file, "r") as f:
            write_numbers = f['scales']['write_number'][()]
            sim_times = f['scales']['sim_time'][()]
        for wn, st in zip(write_numbers, sim_times):
            if cutoff is not None and st > cutoff:
                continue
            if stride > 1 and (int(wn) - 1) % stride != 0:
                continue
            global_t = st + offset
            if start_time is not None and global_t < start_time:
                continue
            if end_time is not None and global_t > end_time:
                continue
            entries.append((int(wn), float(st)))
    entries.sort(key=lambda e: e[1])
    return entries


def compute_offsets(dirs):
    """
    Cumulative time offset for each dir in a restart chain: 0 for the first,
    then the running sum of each preceding run's final-checkpoint sim_time
    (the true sim_time each next run actually restarted from). Used to make
    the burned-in "t = ..." frame titles continuous across a combined movie,
    since each run's own frames are rendered with local, zero-based sim_time.
    """
    offsets = []
    cum = 0.0
    for i, dir in enumerate(dirs):
        offsets.append(cum)
        if i < len(dirs) - 1:
            t = get_checkpoint_final_time(dir)
            cum += t if t is not None else 0.0
    return offsets


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == '--offsets':
        dirs = [d.rstrip('/') + '/' for d in sys.argv[2:]]
        for off in compute_offsets(dirs):
            print(off)
        return

    argv = sys.argv[1:]
    stride = 1
    start_time = None
    end_time = None
    for a in argv:
        if a.startswith('--stride='):
            stride = int(a.split('=', 1)[1])
        elif a.startswith('--start-time='):
            start_time = float(a.split('=', 1)[1])
        elif a.startswith('--end-time='):
            end_time = float(a.split('=', 1)[1])
    argv = [a for a in argv if not (a.startswith('--stride=') or a.startswith('--start-time=') or a.startswith('--end-time='))]

    if len(argv) < 3:
        print("usage: python3 merge_slice_frames.py <plane_tag> <staging_dir> [--stride=N] [--start-time=T] [--end-time=T] <dir1> [<dir2> ...]")
        print("       python3 merge_slice_frames.py --offsets <dir1> <dir2> [<dir3> ...]")
        sys.exit(1)

    plane = argv[0]
    staging_dir = argv[1]
    dirs = [d.rstrip('/') + '/' for d in argv[2:]]

    os.makedirs(staging_dir, exist_ok=True)

    offsets = compute_offsets(dirs)

    frame_paths = []
    for i, dir in enumerate(dirs):
        cutoff = None
        if i < len(dirs) - 1:
            cutoff = get_checkpoint_final_time(dir)
            if cutoff is None:
                print("warning: no checkpoints found in {}, not trimming".format(dir))
        for wn, st in frame_list(dir, cutoff, stride, offsets[i], start_time, end_time):
            frame_path = "{}{}/{}_{:06d}.png".format(dir, plane, plane, wn)
            if os.path.exists(frame_path):
                frame_paths.append(frame_path)
            else:
                print("warning: missing frame {}".format(frame_path))

    if not frame_paths:
        print("no frames found for plane '{}' across: {}".format(plane, dirs))
        sys.exit(1)

    for existing in glob.glob(os.path.join(staging_dir, '*.png')):
        os.remove(existing)

    for n, src in enumerate(frame_paths, start=1):
        dst = os.path.join(staging_dir, "{:06d}.png".format(n))
        os.symlink(os.path.abspath(src), dst)

    print("staged {} frames for plane '{}' in {}".format(len(frame_paths), plane, staging_dir))

if __name__ == "__main__":
    main()
