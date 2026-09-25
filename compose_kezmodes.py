import sys
import os
import glob
import multiprocessing as mp
import subprocess
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from merge_slice_frames import get_checkpoint_final_time, compute_offsets, frame_list

# For each slice-plot frame (as produced by plot_clean.py / staged in time
# order by merge_slice_frames.py), builds a composite image: the frame on the
# left, and the run's kezmodes plot (ke_mode1..10 vs time, spanning the WHOLE
# chain, same trim/offset convention as plot_energies_combined.py and
# merge_slice_frames.py) on the right, with a vertical line marking that
# frame's true simulation time. The base kezmodes plot is rendered once and
# reused -- only the marker line is redrawn per frame -- to keep this cheap
# even over thousands of frames.
#
# Usage:
#   python3 compose_kezmodes.py <plane_tag> <staging_dir> [--stride=N] <dir1> [<dir2> ...]

def color10(label):
    ki = int(label[7:])
    return {
        1: 'black', 2: 'grey', 3: 'lime', 4: 'pink', 5: 'purple',
        6: 'navy', 7: 'green', 8: 'orange', 9: 'yellow', 10: 'brown',
    }[ki]

# be_x/y/z + ke_x/y/z, same six tasks and same colors as
# plot_energies_combined.py's plot_energies(prefix='') -- the source of
# energies__combined.png.
ENERGY_TASKS = ['be_y', 'be_z', 'be_x', 'ke_y', 'ke_z', 'ke_x']

def color_energy(label):
    return {
        'be_y': 'deepskyblue', 'be_z': 'orange', 'be_x': 'green',
        'ke_y': 'magenta', 'ke_z': 'purple', 'ke_x': 'brown',
    }[label]

def label_energy(label):
    return {
        'be_y': r"$0.5\langle b_y^2 \rangle$",
        'be_z': r"$0.5\langle b_z^2 \rangle$",
        'be_x': r"$0.5\langle b_x^2 \rangle$",
        'ke_y': r"$0.5\langle u_y^2 \rangle$",
        'ke_z': r"$0.5\langle u_z^2 \rangle$",
        'ke_x': r"$0.5\langle u_x^2 \rangle$",
    }[label]

def load_series_across_chain(dirs, tasks):
    """`tasks` + sim_time across the chain, trimmed to each run's own final
    checkpoint (except the last) and time-shifted to one continuous series --
    same convention as plot_energies_combined.py's load_series."""
    series = {t: [] for t in tasks}
    series['sim_time'] = []
    offsets = compute_offsets(dirs)
    for i, dir in enumerate(dirs):
        cutoff = get_checkpoint_final_time(dir) if i < len(dirs) - 1 else None
        offset = offsets[i]
        files = glob.glob("{}scalars/*.h5".format(dir))
        files.sort(key=lambda f: int(''.join(filter(str.isdigit, os.path.basename(f)))))
        for file in files:
            with h5py.File(file, "r") as f:
                sim_times = f['scales']['sim_time'][()]
                mask = np.ones(sim_times.shape, dtype=bool)
                if cutoff is not None:
                    mask &= sim_times <= cutoff
                if not mask.any():
                    continue
                series['sim_time'].append(sim_times[mask] + offset)
                for t in tasks:
                    series[t].append(f['tasks'][t][()].ravel()[mask])
    for k in series:
        series[k] = np.concatenate(series[k]) if series[k] else np.array([])
    return series

def load_kezmodes_series(dirs):
    tasks = ['ke_mode{}'.format(k) for k in range(1, 11)]
    return load_series_across_chain(dirs, tasks), tasks

def load_energy_series(dirs):
    return load_series_across_chain(dirs, ENERGY_TASKS), ENERGY_TASKS

def render_base_plot(series, tasks, frame_height_px, aspect=8 / 5, dpi=800):
    # Render at the SAME dpi as plot_clean.py's slice frames (800) and at
    # exactly the frame's pixel height, rather than rendering small and
    # raster-resizing up afterward. Two benefits: no blur from upscaling a
    # low-res render, and -- since matplotlib's default font sizes are used
    # on both sides and dpi now matches -- the kezmodes panel's title/axis/
    # tick text comes out the same point-for-point size as the slice
    # frame's, with no compensating font-size math needed.
    figsize = (frame_height_px / dpi * aspect, frame_height_px / dpi)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    for label in tasks:
        ax.plot(series['sim_time'], series[label], label=label, color=color10(label))
    ax.legend(framealpha=0.0, loc='best', fontsize=7)
    ax.set_ylim(bottom=0)
    ax.set_xlim(0, series['sim_time'][-1])
    ax.set_xlabel("time")
    ax.set_title("ke z-modes")
    fig.tight_layout()
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.asarray(fig.canvas.buffer_rgba()).reshape(h, w, 4)
    img = Image.fromarray(buf, 'RGBA').convert('RGB')
    xlim = ax.get_xlim()
    bbox = ax.get_window_extent()
    plt.close(fig)
    return img, xlim, bbox

def time_to_pixel_x(t, xlim, bbox):
    frac = (t - xlim[0]) / (xlim[1] - xlim[0]) if xlim[1] != xlim[0] else 0.0
    frac = min(max(frac, 0.0), 1.0)
    return bbox.x0 + frac * (bbox.x1 - bbox.x0)

def with_marker(base_img, x_pixel):
    img = base_img.copy()
    draw = ImageDraw.Draw(img)
    draw.line([(x_pixel, 0), (x_pixel, img.height)], fill=(255, 0, 0), width=2)
    return img

def worker_count():
    """Cores to composite with: inside a Slurm allocation, use what it
    granted; otherwise fall back to the node's actual usable core count.
    Deliberately uses `nproc` rather than os.cpu_count()/sched_getaffinity
    for that fallback -- on a shared login node the affinity mask can claim
    e.g. 32 cores while a cgroup CPU-bandwidth quota actually caps the
    process to a fraction of one core; os.cpu_count() ignores that quota and
    forking that many workers just oversubscribes a throttled login node,
    while `nproc` (coreutils, cgroup-aware) reports what's really usable."""
    for var in ('SLURM_NTASKS', 'SLURM_CPUS_ON_NODE'):
        v = os.environ.get(var)
        if v:
            return int(v)
    try:
        return int(subprocess.check_output(['nproc']).strip())
    except Exception:
        return os.cpu_count() or 1

_worker_state = {}

def _init_worker(base_img, xlim, bbox, staging_dir):
    # Runs once per forked worker process; base_img/xlim/bbox are inherited
    # via fork's copy-on-write rather than re-pickled per frame.
    _worker_state['base_img'] = base_img
    _worker_state['xlim'] = xlim
    _worker_state['bbox'] = bbox
    _worker_state['staging_dir'] = staging_dir

def _composite_one(args):
    n, frame_path, t = args
    marked = with_marker(_worker_state['base_img'], time_to_pixel_x(t, _worker_state['xlim'], _worker_state['bbox']))
    composite = compose(frame_path, marked)
    composite.save(os.path.join(_worker_state['staging_dir'], "{:06d}.png".format(n)))
    return n

def compose(frame_path, kezmodes_img):
    # kezmodes_img is already rendered at exactly frame height (see
    # render_base_plot), so this is a plain side-by-side paste -- no resize,
    # which is what kept the marker line's apparent width and the panel's
    # font size drifting/mismatched before.
    frame = Image.open(frame_path).convert('RGB')
    canvas = Image.new('RGB', (frame.width + kezmodes_img.width, frame.height), (255, 255, 255))
    canvas.paste(frame, (0, 0))
    canvas.paste(kezmodes_img, (frame.width, 0))
    return canvas


def main():
    argv = sys.argv[1:]
    stride = 1
    for a in argv:
        if a.startswith('--stride='):
            stride = int(a.split('=', 1)[1])
    argv = [a for a in argv if not a.startswith('--stride=')]

    if len(argv) < 3:
        print("usage: python3 compose_kezmodes.py <plane_tag> <staging_dir> [--stride=N] <dir1> [<dir2> ...]")
        sys.exit(1)

    plane = argv[0]
    staging_dir = argv[1]
    dirs = [d.rstrip('/') + '/' for d in argv[2:]]

    os.makedirs(staging_dir, exist_ok=True)
    for existing in glob.glob(os.path.join(staging_dir, '*.png')):
        os.remove(existing)

    series, tasks = load_kezmodes_series(dirs)
    if len(series['sim_time']) == 0:
        print("no ke_mode data found across: {}".format(dirs))
        sys.exit(1)

    offsets = compute_offsets(dirs)
    frame_entries = []
    for i, dir in enumerate(dirs):
        cutoff = get_checkpoint_final_time(dir) if i < len(dirs) - 1 else None
        offset = offsets[i]
        for wn, st in frame_list(dir, cutoff, stride):
            frame_path = "{}{}/{}_{:06d}.png".format(dir, plane, plane, wn)
            if os.path.exists(frame_path):
                frame_entries.append((frame_path, st + offset))
            else:
                print("warning: missing frame {}".format(frame_path))

    if not frame_entries:
        print("no frames found for plane '{}' across: {}".format(plane, dirs))
        sys.exit(1)

    # Render the kezmodes panel to match this plane's actual frame height
    # (see render_base_plot) -- read it off the first real frame.
    frame_height_px = Image.open(frame_entries[0][0]).height
    base_img, xlim, bbox = render_base_plot(series, tasks, frame_height_px)

    total = len(frame_entries)
    nprocs = min(worker_count(), total)
    print("compositing {} frames for plane '{}' across {} processes...".format(total, plane, nprocs), flush=True)
    progress_every = max(1, total // 20)
    tasks = [(n, frame_path, t) for n, (frame_path, t) in enumerate(frame_entries, start=1)]
    done = 0
    ctx = mp.get_context('fork')
    with ctx.Pool(processes=nprocs, initializer=_init_worker, initargs=(base_img, xlim, bbox, staging_dir)) as pool:
        for _ in pool.imap_unordered(_composite_one, tasks):
            done += 1
            if done % progress_every == 0 or done == total:
                print("  composited {}/{} frames".format(done, total), flush=True)

    print("staged {} kezmodes-composite frames for plane '{}' in {}".format(len(frame_entries), plane, staging_dir), flush=True)

if __name__ == "__main__":
    main()
