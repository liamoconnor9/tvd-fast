import sys
import os
import glob
import multiprocessing as mp
from configparser import ConfigParser

import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
plt.ioff()

from dedalus.extras import plot_tools

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from merge_slice_frames import get_checkpoint_final_time, compute_offsets, natural_sorted_h5
from compose_kezmodes import (
    load_kezmodes_series, load_energy_series, color10, color_energy, label_energy, worker_count,
)

# Renders the slice-plot panel(s) (as plot_clean.py does, via dedalus's
# plot_tools.plot_bot) and a side panel -- either the run's kezmodes plot
# (ke_mode1..10 vs time) or its energies plot (be_x/y/z, ke_x/y/z vs time,
# same six tasks/colors/log-scale as energies__combined.png), spanning the
# whole chain, with a vertical line marking the current frame's time -- into
# ONE matplotlib figure per frame, instead of raster-pasting two separately
# rendered images together (that's what compose_kezmodes.py did, and it
# couldn't make the two panels' plot frames line up, or bound the marker line
# to the side panel's axes, without reverse-engineering pixel offsets).
# Building both panels in the same figure means:
#   - the side panel's axes box (its "frame") is placed at exactly the same
#     figure-fraction y0/y1 as the leftmost slice panel's axes box, so the two
#     frames line up;
#   - the marker is a real matplotlib Line2D (ax.axvline) drawn in the side
#     panel's axes, so it's clipped to that axes' frame for free and can be
#     given the exact same linewidth as the lines it's drawn among;
#   - the "t = ..." suptitle is centered over the whole combined figure.
#
# Usage:
#   python3 plot_kezmodes_combined.py <plane_tag> <staging_dir> [--stride=N] [--dark] [--kezmodes] [--energies] <dir1> [<dir2> ...]
# e.g.
#   python3 plot_kezmodes_combined.py cleany /tmp/stage123 Ro3p5_Lz8pi_Ny16_Rm175 Ro3p5_Lz8pi_Ny16_Rm175_RSTRT1
#   python3 plot_kezmodes_combined.py cleany /tmp/stage123 --energies Ro3p5_Lz8pi_Ny16_Rm175
#   python3 plot_kezmodes_combined.py cleany /tmp/stage123 --kezmodes --energies Ro3p5_Lz8pi_Ny16_Rm175
#
# With both --kezmodes and --energies, the two panels are stacked vertically
# (energies on top, kezmodes below) in the same right-hand column rather than
# each getting its own video: they split the same [frame_y0, frame_y1] span
# the single panel would otherwise occupy (so the column's top/bottom edges
# don't move), share one x-axis (only the bottom/kezmodes panel gets x tick
# labels and the "time" label), and share one xlim so the two line up.

PANEL_LINEWIDTH = 1.2
PANEL_ASPECT = (1 + 5 ** 0.5) / 2  # golden ratio: side panel width : height
PANEL_GAP_FRAC = 0.03  # gap between colormesh panels and the side panel, as a fraction of the frame height
STACK_GAP_FRAC = 0.05  # vertical gap between stacked energies/kezmodes panels, as a fraction of the frame height
DPI = 800  # match plot_clean.py's slice-frame dpi


# ---- small helpers duplicated from plot_clean.py, so this script has no
# import-time coupling to its MPI/global-state (this runs as a plain
# multiprocessing pool, not under mpirun). ----

def dark_colormaps():
    from cmap import Colormap
    return Colormap('vanimo').to_mpl(), Colormap('berlin').to_mpl()

def place_offset_label(cb_ax, oom):
    """Put the 'x10^oom' scale factor just to the right of this colorbar's
    own title (whatever it is -- short like "b_y" or long like the decoy
    label), instead of matplotlib's default top-right placement (which
    collides with the title) or a hardcoded offset (which wouldn't
    generalize across title lengths). Requires a draw so the title's
    rendered extent is known."""
    from matplotlib.transforms import offset_copy
    fig = cb_ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    title_text = cb_ax.xaxis.label
    bbox = title_text.get_window_extent(renderer=renderer)
    inv = cb_ax.transAxes.inverted()
    x, y = inv.transform((bbox.x1, (bbox.y0 + bbox.y1) / 2))
    label_transform = offset_copy(cb_ax.transAxes, fig=fig, x=6, y=0, units='points')
    fontsize = cb_ax.get_xticklabels()[0].get_fontsize()
    cb_ax.text(x, y, r"$\times 10^{{{}}}$".format(oom),
               transform=label_transform, ha='left', va='center', fontsize=fontsize)

def set_centered_ticks(cb_ax):
    """3 ticks (min, 0, max) at 2 significant figures instead of dedalus's
    default 5 at 3 sig figs -- with even_scale=True and no explicit clim,
    the auto 5-tick colorbar labels overlap and collide on these narrow
    subplots. clim is already symmetric about 0 (that's what even_scale
    does), so this just relabels it clearly without changing the
    (already zero-centered) color scale itself."""
    lo, hi = cb_ax.get_xlim()
    cb_ax.set_xticks([lo, 0, hi])
    oom = int(np.floor(np.log10(hi))) if hi > 0 else 0
    scale = 10.0 ** oom

    def fmt(v):
        return "0" if v == 0 else "{:.0f}".format(v / scale)

    # 1 sig fig -- short enough that these fit at the plot's regular tick
    # font size without the shrink-to-fit this used to need at 2 sig figs.
    cb_ax.set_xticklabels([fmt(lo), fmt(0), fmt(hi)])
    cb_ax.xaxis.get_offset_text().set_visible(False)
    if oom != 0:
        place_offset_label(cb_ax, oom)

def format_z_axis_pi(paxes, Lz):
    n_pi = int(round(Lz / np.pi))
    if n_pi < 1:
        return
    step_k = max(1, round(n_pi / 5))
    ticks_k = sorted(set(range(0, n_pi, step_k)) | {n_pi})
    tick_positions = [k * np.pi for k in ticks_k]
    paxes.set_ylim(top=Lz)

    def pi_formatter(z, pos):
        k = z / np.pi
        k_int = int(round(k))
        if abs(k - k_int) > 1e-6:
            return ""
        if k_int == 0:
            return "0"
        if k_int == 1:
            return r"$\pi$"
        if k_int == -1:
            return r"$-\pi$"
        return r"${}\pi$".format(k_int)

    paxes.yaxis.set_major_locator(ticker.FixedLocator(tick_positions))
    paxes.yaxis.set_major_formatter(ticker.FuncFormatter(pi_formatter))

def hide_shared_yaxis(paxes):
    paxes.set_ylabel('')
    paxes.tick_params(axis='y', labelleft=False)


def read_config(dir):
    config = ConfigParser()
    config.read("{}options.cfg".format(dir))
    Ly = eval(config.get('parameters', 'Ly'))
    Lz = eval(config.get('parameters', 'Lz'))
    Lx = eval(config.get('parameters', 'Lx'))
    isHydro = config.getboolean('parameters', 'isHydro')
    is2d = config.getboolean('parameters', 'is2d')
    return Ly, Lz, Lx, isHydro, is2d


def frame_entries(dir, cutoff, stride, offset=0.0, start_time=None, end_time=None):
    """[(h5_path, h5_index, write_number, sim_time), ...] for slicepoints
    writes in `dir`, chronologically ordered, trimmed to sim_time <= cutoff,
    and subsampled to every `stride`-th write -- same predicate as
    merge_slice_frames.frame_list, but keeping (file, index) so the raw
    pcolor data can be read directly instead of a pre-rendered PNG path.

    start_time/end_time bound the combined/global time (st + offset), same
    meaning as in merge_slice_frames.frame_list / plot_clean.py."""
    entries = []
    for file in natural_sorted_h5("{}slicepoints/*.h5".format(dir)):
        with h5py.File(file, "r") as f:
            write_numbers = f['scales']['write_number'][()]
            sim_times = f['scales']['sim_time'][()]
        for idx, (wn, st) in enumerate(zip(write_numbers, sim_times)):
            if cutoff is not None and st > cutoff:
                continue
            if stride > 1 and (int(wn) - 1) % stride != 0:
                continue
            global_t = st + offset
            if start_time is not None and global_t < start_time:
                continue
            if end_time is not None and global_t > end_time:
                continue
            entries.append((file, idx, int(wn), float(st)))
    entries.sort(key=lambda e: e[3])
    return entries


def plane_geometry(normal_dir, ary, arz):
    if normal_dir == 'x':
        image = plot_tools.Box(2, 2 * arz / ary)
        image_axes = (1, 2)
        data_slices_tail = (slice(None), slice(None), 0)
    elif normal_dir == 'y':
        image = plot_tools.Box(2 / arz, 2)
        image_axes = (3, 2)
        data_slices_tail = (0, slice(None), slice(None))
    else:
        raise ValueError("only normal_dir 'x'/'y' are used by --kezmodes")
    return image, image_axes, data_slices_tail


def draw_panel(panel_ax, panel_mode, panel_series, panel_tasks, dark, show_xlabel=True, xlim_max=None):
    if panel_mode == 'energies':
        for label in panel_tasks:
            linestyle = 'dashed' if label.startswith('ke') else 'solid'
            panel_ax.plot(panel_series['sim_time'], panel_series[label], label=label_energy(label),
                           color=color_energy(label), linewidth=PANEL_LINEWIDTH, linestyle=linestyle)
        panel_ax.legend(framealpha=0.0, loc='best', fontsize=7)
        panel_ax.set_yscale('log')
        panel_ax.set_ylabel('energy')
    else:
        for label in panel_tasks:
            panel_ax.plot(panel_series['sim_time'], panel_series[label], label=label,
                           color=color10(label), linewidth=PANEL_LINEWIDTH)
        panel_ax.legend(framealpha=0.0, loc='best', fontsize=7)
        panel_ax.set_ylim(bottom=0)
    panel_ax.set_xlim(0, xlim_max if xlim_max is not None else panel_series['sim_time'][-1])
    if show_xlabel:
        panel_ax.set_xlabel("time")
    else:
        panel_ax.tick_params(axis='x', labelbottom=False)
    panel_ax.yaxis.tick_right()
    panel_ax.yaxis.set_label_position('right')


def style_dark_panel(panel_ax):
    panel_ax.tick_params(colors='white')
    for spine in panel_ax.spines.values():
        spine.set_color('white')
    panel_ax.xaxis.label.set_color('white')
    panel_ax.yaxis.label.set_color('white')


def build_composite_frame(cfg, dark, panel_mode, panels, file, h5_index, sim_time, decoy):
    normal_dir = cfg['normal_dir']
    isHydro = cfg['isHydro']
    Lz = cfg['Lz']
    cmap_name = cfg['cmap_name']
    cmap_berlin = cfg['cmap_berlin']
    image, image_axes, data_slices_tail = plane_geometry(normal_dir, cfg['ary'], cfg['arz'])

    if not isHydro:
        nrows, ncols = 1, 2
        tasks = ['vy', 'by']
    else:
        nrows, ncols = 1, 1
        tasks = ['vy']
    tasks = [t + '_mid' + normal_dir for t in tasks]
    if normal_dir == 'y' and decoy:
        nrows, ncols = 1, 3
        tasks += ['decoy.ey']

    pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
    margin = plot_tools.Frame(0.3, 0.2, 0.1, 0.1)
    scale = 2.0

    mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
    fig = mfig.figure
    vmin_u, vmax_u = -1, 1

    data_slices = (h5_index,) + data_slices_tail
    frame_y0 = frame_y1 = None
    title_ref_ax = None
    with h5py.File(file, 'r') as f:
        for n, task in enumerate(tasks):
            i, j = divmod(n, ncols)
            axes = mfig.add_axes(i, j, [0, 0, 1, 1])
            dset = f['tasks'][task]

            if 'b' in task:
                pl = plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes,
                                          title=r"$b_y$", even_scale=True, cmap=cmap_name)
                cb_ax = pl[1]
                set_centered_ticks(cb_ax)
            elif 'decoy' in task:
                pl = plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes,
                                          title=r"$((\mathbf{B_0}\cdot\nabla) \mathbf{B_0})\cdot \mathbf{\hat{e_y}}$",
                                          even_scale=True, cmap=cmap_name)
                cb_ax = pl[1]
                set_centered_ticks(cb_ax)
            else:
                if dark:
                    if normal_dir == 'y':
                        pl = plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes,
                                                  title=r"$u_y$", even_scale=True,
                                                  clim=(vmin_u, vmax_u), cmap=cmap_berlin)
                        cb_ax = pl[1]
                        cb_ax.set_xticks([vmin_u, 0, vmax_u])
                    else:
                        pl = plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes,
                                                  title=r"$u_y$", even_scale=True, cmap=cmap_berlin)
                        cb_ax = pl[1]
                        set_centered_ticks(cb_ax)
                else:
                    pl = plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes,
                                              title=r"$u_y$", even_scale=True, clim=(vmin_u, vmax_u))
                    cb_ax = pl[1]
                    cb_ax.set_xticks([vmin_u, 0, vmax_u])

            paxes = pl[0]
            if n == 0:
                format_z_axis_pi(paxes, Lz)
                pos = paxes.get_position()
                frame_y0, frame_y1 = pos.y0, pos.y1
                # Reference point for the "t = ..." suptitle: put it at the
                # same height as this column's title label (e.g. "$u_y$"),
                # not up in the empty margin strip above it.
                title_ref_ax = cb_ax
            else:
                hide_shared_yaxis(paxes)

            if dark:
                cb_ax.tick_params(colors='white')
                for label in cb_ax.get_xticklabels() + cb_ax.get_yticklabels():
                    label.set_color('white')

    # --- widen the canvas and append the side panel, frame-aligned to the
    # leftmost slice panel's own axes box (frame_y0/frame_y1). Existing axes
    # are re-expressed as fractions of the wider canvas (x0/width scaled, y
    # untouched) so their absolute on-page size/position doesn't change --
    # this just frees up real estate on the right. ---
    old_w, old_h = fig.get_size_inches()
    frame_height_in = (frame_y1 - frame_y0) * old_h
    gap_in = PANEL_GAP_FRAC * frame_height_in
    right_pad_in = 0.12 * frame_height_in
    panel_width_in = frame_height_in * PANEL_ASPECT
    new_w = old_w + gap_in + panel_width_in + right_pad_in

    scale_x = old_w / new_w
    for ax in fig.axes:
        p = ax.get_position()
        ax.set_position([p.x0 * scale_x, p.y0, p.width * scale_x, p.height])
    fig.set_size_inches(new_w, old_h, forward=True)

    panel_x0 = (old_w + gap_in) / new_w
    panel_w = panel_width_in / new_w

    if panel_mode == 'both':
        kez_series, kez_tasks = panels['kezmodes']
        energy_series, energy_tasks = panels['energies']
        end_time = max(kez_series['sim_time'][-1], energy_series['sim_time'][-1])
        frame_h = frame_y1 - frame_y0
        stack_gap = STACK_GAP_FRAC * frame_h
        sub_h = (frame_h - stack_gap) / 2
        # kezmodes on the bottom (its own bottom edge stays at frame_y0),
        # energies on top (its own top edge stays at frame_y1) -- together
        # they span exactly the same [frame_y0, frame_y1] the single panel
        # used to occupy.
        kez_ax = fig.add_axes([panel_x0, frame_y0, panel_w, sub_h])
        energy_ax = fig.add_axes([panel_x0, frame_y0 + sub_h + stack_gap, panel_w, sub_h])

        draw_panel(energy_ax, 'energies', energy_series, energy_tasks, dark,
                   show_xlabel=False, xlim_max=end_time)
        draw_panel(kez_ax, 'kezmodes', kez_series, kez_tasks, dark,
                   show_xlabel=True, xlim_max=end_time)

        for panel_ax in (kez_ax, energy_ax):
            # Bounded to this axes' own frame (clip_on defaults True) and
            # exactly as wide as the lines it's drawn among.
            panel_ax.axvline(sim_time, color='red', linewidth=PANEL_LINEWIDTH, zorder=5)
            if dark:
                style_dark_panel(panel_ax)
    else:
        panel_series, panel_tasks = panels[panel_mode]
        panel_ax = fig.add_axes([panel_x0, frame_y0, panel_w, frame_y1 - frame_y0])

        draw_panel(panel_ax, panel_mode, panel_series, panel_tasks, dark)
        panel_ax.axvline(sim_time, color='red', linewidth=PANEL_LINEWIDTH, zorder=5)
        if dark:
            style_dark_panel(panel_ax)

    # Same height as the colorbar title labels (e.g. "$u_y$"/"$b_y$"), not
    # the empty margin strip above them -- read the reference label's actual
    # rendered position rather than guessing at a fraction of the margin.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_bbox = title_ref_ax.xaxis.label.get_window_extent(renderer=renderer)
    title_y = fig.transFigure.inverted().transform((0, (label_bbox.y0 + label_bbox.y1) / 2))[1]

    title = 't = {:.3f}'.format(sim_time)
    if dark:
        fig.suptitle(title, x=0.5, y=title_y, ha='center', va='center', color='white')
    else:
        fig.suptitle(title, x=0.5, y=title_y, ha='center', va='center')

    return fig


_worker_state = {}

def _init_worker(cfg, dark, panel_mode, panels, staging_dir):
    if dark:
        plt.style.use('dark_background')
    _worker_state['cfg'] = cfg
    _worker_state['dark'] = dark
    _worker_state['panel_mode'] = panel_mode
    _worker_state['panels'] = panels
    _worker_state['staging_dir'] = staging_dir

def _composite_one(args):
    n, file, h5_index, sim_time, decoy = args
    fig = build_composite_frame(_worker_state['cfg'], _worker_state['dark'],
                                 _worker_state['panel_mode'], _worker_state['panels'],
                                 file, h5_index, sim_time, decoy)
    fig.savefig(os.path.join(_worker_state['staging_dir'], "{:06d}.png".format(n)),
                dpi=DPI, facecolor=fig.get_facecolor())
    plt.close(fig)
    return n


def main():
    argv = sys.argv[1:]
    stride = 1
    dark = False
    do_energies = False
    do_kezmodes = False
    start_time = None
    end_time = None
    rest = []
    for a in argv:
        if a.startswith('--stride='):
            stride = int(a.split('=', 1)[1])
        elif a.startswith('--start-time='):
            start_time = float(a.split('=', 1)[1])
        elif a.startswith('--end-time='):
            end_time = float(a.split('=', 1)[1])
        elif a in ('--dark', '-d'):
            dark = True
        elif a in ('--energies', '-e'):
            do_energies = True
        elif a in ('--kezmodes', '-k'):
            do_kezmodes = True
        else:
            rest.append(a)

    if not do_energies and not do_kezmodes:
        do_kezmodes = True

    if len(rest) < 3:
        print("usage: python3 plot_kezmodes_combined.py <plane_tag> <staging_dir> [--stride=N] [--start-time=T] [--end-time=T] [--dark] [--kezmodes] [--energies] <dir1> [<dir2> ...]")
        sys.exit(1)

    plane = rest[0]
    normal_dir = {'cleany': 'y', 'cleanx': 'x'}.get(plane)
    if normal_dir is None:
        print("unknown plane_tag '{}' (expected cleany or cleanx)".format(plane))
        sys.exit(1)
    staging_dir = rest[1]
    dirs = [d.rstrip('/') + '/' for d in rest[2:]]

    os.makedirs(staging_dir, exist_ok=True)
    for existing in glob.glob(os.path.join(staging_dir, '*.png')):
        os.remove(existing)

    Ly, Lz, Lx, isHydro, is2d = read_config(dirs[0])
    if normal_dir == 'x' and is2d:
        print("2D run: no cleanx frames to render")
        sys.exit(1)

    if dark:
        cmap_name, cmap_berlin = dark_colormaps()
    else:
        cmap_name, cmap_berlin = "PiYG", None

    cfg = {
        'normal_dir': normal_dir,
        'isHydro': isHydro,
        'Lz': Lz,
        'ary': Ly / Lx,
        'arz': Lz / Lx,
        'cmap_name': cmap_name,
        'cmap_berlin': cmap_berlin,
    }

    panel_mode = 'both' if (do_kezmodes and do_energies) else ('energies' if do_energies else 'kezmodes')

    panels = {}
    if do_kezmodes:
        kez_series, kez_tasks = load_kezmodes_series(dirs)
        if len(kez_series['sim_time']) == 0:
            print("no ke_mode data found across: {}".format(dirs))
            sys.exit(1)
        panels['kezmodes'] = (kez_series, kez_tasks)
    if do_energies:
        energy_series, energy_tasks = load_energy_series(dirs)
        if len(energy_series['sim_time']) == 0:
            print("no energy data found across: {}".format(dirs))
            sys.exit(1)
        panels['energies'] = (energy_series, energy_tasks)

    offsets = compute_offsets(dirs)
    frame_tasks = []
    n = 0
    for i, dir in enumerate(dirs):
        cutoff = get_checkpoint_final_time(dir) if i < len(dirs) - 1 else None
        offset = offsets[i]
        decoy = 'decoy' in os.path.basename(dir.rstrip('/'))
        for file, h5_index, wn, st in frame_entries(dir, cutoff, stride, offset, start_time, end_time):
            n += 1
            frame_tasks.append((n, file, h5_index, st + offset, decoy))

    if not frame_tasks:
        print("no frames found for plane '{}' across: {}".format(plane, dirs))
        sys.exit(1)

    total = len(frame_tasks)
    nprocs = min(worker_count(), total)
    print("rendering {} combined frames for plane '{}' across {} processes...".format(total, plane, nprocs), flush=True)
    progress_every = max(1, total // 20)
    done = 0
    ctx = mp.get_context('fork')
    with ctx.Pool(processes=nprocs, initializer=_init_worker,
                  initargs=(cfg, dark, panel_mode, panels, staging_dir)) as pool:
        for _ in pool.imap_unordered(_composite_one, frame_tasks):
            done += 1
            if done % progress_every == 0 or done == total:
                print("  rendered {}/{} frames".format(done, total), flush=True)

    print("staged {} combined frames for plane '{}' in {}".format(len(frame_tasks), plane, staging_dir), flush=True)

if __name__ == "__main__":
    main()
