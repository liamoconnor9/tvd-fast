from glob import glob
from configparser import ConfigParser
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.transforms import offset_copy
plt.ioff()

from dedalus.extras import plot_tools
import logging
import sys
logger = logging.getLogger(__name__)
import os
path = os.path.dirname(os.path.abspath(__file__))
from mpi4py import MPI
CW = MPI.COMM_WORLD

# Dark-mode colormaps (only built/used when --dark is passed).
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
    fig = cb_ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    title_text = cb_ax.xaxis.label
    bbox = title_text.get_window_extent(renderer=renderer)
    inv = cb_ax.transAxes.inverted()
    x, y = inv.transform((bbox.x1, (bbox.y0 + bbox.y1) / 2))
    label_transform = offset_copy(cb_ax.transAxes, fig=fig, x=6, y=0, units='points')
    cb_ax.text(x, y, r"$\times 10^{{{}}}$".format(oom),
               transform=label_transform, ha='left', va='center', fontsize=6)

def shrink_ticklabels_to_fit(cb_ax, min_fontsize=5):
    """Shrink the tick labels' font size until adjacent ones no longer
    overlap. How much room 3 tick labels have depends on the colorbar's
    own width, which shrinks as Lz grows (the multi-panel box gets
    proportionally narrower/taller) -- so a fixed font size that's fine
    for a short domain can still collide on a long one."""
    fig = cb_ax.figure
    fontsize = cb_ax.get_xticklabels()[0].get_fontsize()
    while fontsize > min_fontsize:
        fig.canvas.draw()
        boxes = [t.get_window_extent(renderer=fig.canvas.get_renderer())
                 for t in cb_ax.get_xticklabels()]
        if all(boxes[i].x1 <= boxes[i + 1].x0 for i in range(len(boxes) - 1)):
            break
        fontsize -= 1
        for t in cb_ax.get_xticklabels():
            t.set_fontsize(fontsize)

def set_centered_ticks(cb_ax):
    """3 ticks (min, 0, max) at 2 significant figures instead of dedalus's
    default 5 at 3 sig figs -- with even_scale=True and no explicit clim,
    the auto 5-tick colorbar labels overlap and collide on these narrow
    subplots. clim is already symmetric about 0 (that's what even_scale
    does), so this just relabels it clearly without changing the
    (already zero-centered) color scale itself. When the data needs a
    order-of-magnitude offset (e.g. values ~1e-5), it's placed beside the
    title via place_offset_label instead of overlapping it; when the data
    is already O(1) (oom == 0, e.g. u_y ~ 1), no offset text is shown at
    all."""
    lo, hi = cb_ax.get_xlim()
    cb_ax.set_xticks([lo, 0, hi])
    oom = int(np.floor(np.log10(hi))) if hi > 0 else 0
    scale = 10.0 ** oom

    def fmt(v):
        return "0.0" if v == 0 else "{:.1f}".format(v / scale)

    cb_ax.set_xticklabels([fmt(lo), fmt(0), fmt(hi)])
    cb_ax.xaxis.get_offset_text().set_visible(False)
    shrink_ticklabels_to_fit(cb_ax)
    if oom != 0:
        place_offset_label(cb_ax, oom)

def format_z_axis_pi(paxes):
    """z spans [0, Lz] with Lz = n*pi for this project's runs -- label the
    shared z-axis (left subplot only; see hide_shared_yaxis) in multiples
    of pi instead of raw radians, always including the top (Lz itself)
    regardless of step size or n. The z grid is periodic (samples
    [0, Lz - dz]), so the axis view never actually reaches Lz on its own;
    extend it to Lz so the top tick has somewhere to be drawn."""
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
    """The right subplot(s) share the same z-axis as the leftmost one --
    drop the repeated label and tick numbers."""
    paxes.set_ylabel('')
    paxes.tick_params(axis='y', labelleft=False)


def plot_plane(filename, start, count, output, normal_dir, tag):
    """Save plot of specified tasks for given range of analysis writes."""

    if normal_dir == 'x':
        image = plot_tools.Box(2, 2 * arz / ary)
        image_axes = (1, 2)
        data_slices_tail = (slice(None), slice(None), 0)

    if normal_dir == 'y':
        image = plot_tools.Box(2 / arz, 2)
        image_axes = (3, 2)
        data_slices_tail = (0, slice(None), slice(None))

    if normal_dir == 'z':
        image = plot_tools.Box(2, 2 / ary)
        image_axes = (1, 3)
        data_slices_tail = (slice(None), 0, slice(None))

    # Plot settings
    scale = 2.0
    dpi = 800
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: '{}_{:06}.png'.format(tag, write)

    if not isHydro:
        nrows, ncols = 1, 2
        tasks = ['vy', 'by']
    else:
        nrows, ncols = 1, 1
        tasks = ['vy']
    tasks = [task + '_mid' + normal_dir for task in tasks]
    if normal_dir == 'y' and 'decoy' in suffix:
        if CW.rank == 0:
            print('up')
        nrows, ncols = 1, 3
        tasks += ['decoy.ey']

    pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
    margin = plot_tools.Frame(0.3, 0.2, 0.1, 0.1)

    # Create multifigure
    mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
    fig = mfig.figure
    vmin_b=-3
    vmax_b=3
    vmin_u=-1
    vmax_u=1

    # Plot writes
    with h5py.File(filename, mode='r') as file:
        for index in range(start, start+count):
            write_number = int(file['scales/write_number'][index])
            if STRIDE > 1 and (write_number - 1) % STRIDE != 0:
                continue
            sim_time = float(file['scales/sim_time'][index]) + TIME_OFFSET
            if START_TIME is not None and sim_time < START_TIME:
                continue
            if END_TIME is not None and sim_time > END_TIME:
                continue
            data_slices = (index, ) + data_slices_tail
            for n, task in enumerate(tasks):
                # Build subfigure axes
                i, j = divmod(n, ncols)
                axes = mfig.add_axes(i, j, [0, 0, 1, 1])
                dset = file['tasks'][task]

                if 'b' in task:
                    pl = plot_tools.plot_bot(
                        dset, image_axes, data_slices,
                        axes=axes, title=r"$b_y$", even_scale=True, cmap=cmap_name
                    )
                    cb_ax = pl[1]
                    set_centered_ticks(cb_ax)
                elif 'decoy' in task:
                    pl = plot_tools.plot_bot(
                        dset, image_axes, data_slices,
                        axes=axes, title=r"$((\mathbf{B_0}\cdot\nabla) \mathbf{B_0})\cdot \mathbf{\hat{e_y}}$",
                        even_scale=True, cmap=cmap_name
                    )
                    cb_ax = pl[1]
                    set_centered_ticks(cb_ax)
                else:
                    if DARK:
                        if normal_dir == 'y':
                            pl = plot_tools.plot_bot(
                                dset, image_axes, data_slices,
                                axes=axes, title=r"$u_y$", even_scale=True, clim=(vmin_u, vmax_u), cmap=cmap_berlin
                            )
                            cb_ax = pl[1]
                            cb_ax.set_xticks([vmin_u, 0, vmax_u])
                        else:
                            pl = plot_tools.plot_bot(
                                dset, image_axes, data_slices,
                                axes=axes, title=r"$u_y$", even_scale=True, cmap=cmap_berlin
                            )
                            cb_ax = pl[1]
                            set_centered_ticks(cb_ax)
                    else:
                        pl = plot_tools.plot_bot(
                            dset, image_axes, data_slices,
                            axes=axes, title=r"$u_y$", even_scale=True, clim=(vmin_u, vmax_u)
                        )
                        cb_ax = pl[1]
                        cb_ax.set_xticks([vmin_u, 0, vmax_u])

                paxes = pl[0]
                if n == 0:
                    format_z_axis_pi(paxes)
                else:
                    hide_shared_yaxis(paxes)

                if DARK:
                    # Force ticks and labels to white
                    cb_ax.tick_params(colors='white')
                    for label in cb_ax.get_xticklabels() + cb_ax.get_yticklabels():
                        label.set_color('white')

            # Add time title
            title = title_func(sim_time)
            title_height = 1 - 0.5 * mfig.margin.top / mfig.fig.y
            if DARK:
                fig.suptitle(title, x=0.42, y=title_height, ha='left', color='white')
            else:
                fig.suptitle(title, x=0.42, y=title_height, ha='left')

            # Save figure
            savename = savename_func(write_number)
            savepath = output.joinpath(savename)
            if (index % 1 == 0):
                if DARK:
                    fig.savefig(str(savepath), dpi=dpi, facecolor=fig.get_facecolor())
                else:
                    fig.savefig(str(savepath), dpi=dpi)
            fig.clear()
    plt.close(fig)


def plot_all(filename, start, count):
    plot_plane(filename, start, count, output_mid_path_cleany, 'y', 'cleany')
    if not is2d:
        plot_plane(filename, start, count, output_mid_path_cleanx, 'x', 'cleanx')
    # plot_plane(filename, start, count, output_mid_path_cleanz, 'z', 'cleanz')


if __name__ == "__main__":

    import pathlib
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync
    global suffix, ar, ary, arz, last_index, isHydro, solver_name, is2d, DARK, cmap_name, cmap_berlin, TIME_OFFSET, STRIDE, START_TIME, END_TIME

    argv = sys.argv[1:]
    DARK = '--dark' in argv or '-d' in argv
    TIME_OFFSET = 0.0
    STRIDE = 1
    START_TIME = None
    END_TIME = None
    for a in argv:
        if a.startswith('--time-offset='):
            TIME_OFFSET = float(a.split('=', 1)[1])
        elif a.startswith('--stride='):
            STRIDE = int(a.split('=', 1)[1])
        elif a.startswith('--start-time='):
            START_TIME = float(a.split('=', 1)[1])
        elif a.startswith('--end-time='):
            END_TIME = float(a.split('=', 1)[1])
    positional = [a for a in argv if not a.startswith('-')]

    if len(positional) > 0:
        suffix = positional[0]
        if suffix[-1] == '/':
            suffix = suffix[:-1]
    else:
        raise

    if DARK:
        plt.style.use('dark_background')
        cmap_name, cmap_berlin = dark_colormaps()
    else:
        cmap_name = "PiYG"
        cmap_berlin = None

    filename = "{}/{}/options.cfg".format(path, suffix)
    config = ConfigParser()
    config.read(str(filename))

    Ly = eval(config.get('parameters','Ly'))
    Lz = eval(config.get('parameters','Lz'))
    Lx = eval(config.get('parameters','Lx'))
    isHydro = config.getboolean('parameters','isHydro')
    is2d = config.getboolean('parameters', 'is2d')
    solver_name = eval(config.get('parameters','SOLVER'))
    if CW.rank == 0:
        print(solver_name)
        print(solver_name == 'ql_fly.py')

    ary = Ly / Lx
    arz = Lz / Lx

    slicepoints = glob("{}/{}/slicepoints/*.h5".format(path, suffix))
    last_index = len(slicepoints)

    # Create output directory if needed
    output_mid_path_cleany=pathlib.Path('{}/{}/cleany'.format(path, suffix))
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_mid_path_cleany.exists():
                output_mid_path_cleany.mkdir()

    output_mid_path_cleanx=pathlib.Path('{}/{}/cleanx'.format(path, suffix))
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_mid_path_cleanx.exists():
                output_mid_path_cleanx.mkdir()

    output_mid_path_cleanz=pathlib.Path('{}/{}/cleanz'.format(path, suffix))
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_mid_path_cleanz.exists():
                output_mid_path_cleanz.mkdir()

    post.visit_writes(slicepoints, plot_all)
