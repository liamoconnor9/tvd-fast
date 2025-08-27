from glob import glob
from docopt import docopt
from configparser import ConfigParser
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
plt.ioff()
from dedalus.extras import plot_tools
import logging
import sys
logger = logging.getLogger(__name__)
import os
path = os.path.dirname(os.path.abspath(__file__))

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
    scale = 2.5
    dpi = 100
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: '{}_{:06}.png'.format(tag, write)
    if not isHydro:
        nrows, ncols = 1, 6
        tasks = ['vy', 'vz', 'vx', 'by', 'bz', 'bx']
        # tasks = ['vy', 'vz', 'vx', 'by', 'bz', 'bx', 'jy', 'jz', 'jx']
    else:
        nrows, ncols = 1, 3
        tasks = ['vy', 'vz', 'vx']
    tasks = [task + '_' + tag + normal_dir for task in tasks]
    if False:
        ncols = 7
        if normal_dir == 'y':
            tasks += tasks + ['bhelicity_midy_midy']
        if normal_dir == 'z':
            tasks += tasks + ['bhelicity_midz_midz']
        if normal_dir == 'y':
            tasks += tasks + ['bhelicity_midz_midz']

    pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
    margin = plot_tools.Frame(0.3, 0.2, 0.1, 0.1)

    # Create multifigure
    mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
    fig = mfig.figure
    # Plot writes
    with h5py.File(filename, mode='r') as file:
        for index in range(start, start+count):
            data_slices = (index, ) + data_slices_tail
            for n, task in enumerate(tasks):
                # Build subfigure axes
                i, j = divmod(n, ncols)
                axes = mfig.add_axes(i, j, [0, 0, 1, 1])
                dset = file['tasks'][task]
                if 'b' in task:
                    plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes, title=task, even_scale=True, cmap="PiYG")
                else:
                    plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes, title=task, even_scale=True)

            # Add time title
            title = title_func(file['scales/sim_time'][index])
            title_height = 1 - 0.5 * mfig.margin.top / mfig.fig.y
            fig.suptitle(title, x=0.48, y=title_height, ha='left')
            # Save figure
            savename = savename_func(file['scales/write_number'][index])
            savepath = output.joinpath(savename)
            if (index % 1 == 0):
                fig.savefig(str(savepath), dpi=dpi)
            fig.clear()
    plt.close(fig)

def plot_all(filename, start, count):
    # if "_s{}.h5".format(last_index) in filename:
        # print('success')
        # return
    plot_plane(filename, start, count, output_mid_path_yz, 'x', 'mid')
    plot_plane(filename, start, count, output_mid_path_zx, 'y', 'mid')
    plot_plane(filename, start, count, output_mid_path_xy, 'z', 'mid')
    # plot_plane(filename, start, count, output_avg_path_yz, 'x', 'avg')
    # plot_plane(filename, start, count, output_avg_path_zx, 'y', 'avg')
    # plot_plane(filename, start, count, output_avg_path_xy, 'z', 'avg')

def yzmean(filename, start, count, output):
    """Save plot of specified tasks for given range of analysis writes."""
    if "_s{}.h5".format(last_index) in filename:
        print('success')
        # return

    # Plot settings
    normal_dir = 'z'
    scale = 2.5
    dpi = 100
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: 'mid_{:06}.png'.format(write)
    # Layout
    # if (round(ary) > 1):
    if not isHydro:
        nrows, ncols = 2, 3
        tasks = ['vy_avg', 'by_avg', 'jy_avg', 'vz_avg', 'bz_avg', 'jz_avg']
    else:
        nrows, ncols = 2, 1
        tasks = ['vy_avg', 'vz_avg']


    # Plot writes
    with h5py.File(filename, mode='r') as file:
        # print()
        for key in file['scales'].keys():
            if 'x_hash' in key:
                x = file['scales'][key][()]

        for index in range(start, start+count):
            fig, axes = plt.subplots(nrows, ncols, sharex=True, sharey=False, layout='constrained')
            for n, task in enumerate(tasks):
                # Build subfigure axes
                i, j = divmod(n, ncols)
                # Call plotting helper (dset axes: [t, x, y, z])
                dset_vec = file['tasks'][task][()][index, 0, 0, :]
                # image_axes = (1, 3)
                # data_slices = (index, slice(None), 0, slice(None))
                # plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes, title=task, even_scale=True)
                axes[i, j].plot(x, dset_vec, linewidth=3, color='purple')
                axes[i, j].set_title(task)
                if i == 1:
                    axes[i, j].set_xlabel('x')
                axes[i, j].set_xlim(-0.5, 0.5)

            # # Add time title
            sim_time = file['scales/sim_time'][index]
            fig.suptitle('t = {:.3f}'.format(sim_time))
            # # Save figure
            savename = savename_func(file['scales/write_number'][index])
            savepath = output.joinpath(savename)
            if (index % 1 == 0):
                plt.savefig(str(savepath), dpi=dpi)
            plt.close()


if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    if len(sys.argv) > 1:
        suffix = sys.argv[1]
        if suffix[-1] == '/':
            suffix = suffix[:-1]
    else:
        raise
    # sys.exit()
    
    # args = docopt(__doc__)

    # output_path = pathlib.Path(args['--output']).absolute()
    # dir = args['--dir']
    # suffix = args['--suffix']
    filename = "{}/{}/options.cfg".format(path, suffix)
    config = ConfigParser()
    config.read(str(filename))

    global ar, ary, arz, last_index, isHydro
    Ly = eval(config.get('parameters','Ly'))
    Lz = eval(config.get('parameters','Lz'))
    Lx = eval(config.get('parameters','Lx'))
    isHydro = config.getboolean('parameters','isHydro')

    ary = Ly / Lx
    arz = Lz / Lx 

    data_dirs = glob("{}/{}//".format(path, suffix))
    for data_dir in data_dirs:
        slicepoints = glob("{}/{}/slicepoints/*.h5".format(path, suffix))
        last_index = len(slicepoints)

        # Create output directory if needed
        output_mid_path_yz=pathlib.Path('{}/{}/mid_yz'.format(path, suffix))
        output_mid_path_zx=pathlib.Path('{}/{}/mid_zx'.format(path, suffix))
        output_mid_path_xy=pathlib.Path('{}/{}/mid_xy'.format(path, suffix))
        # output_avg_path_yz=pathlib.Path('{}/{}/avg_yz'.format(path, suffix))
        # output_avg_path_zx=pathlib.Path('{}/{}/avg_zx'.format(path, suffix))
        # output_avg_path_xy=pathlib.Path('{}/{}/avg_xy'.format(path, suffix))
        # output_path_avg=pathlib.Path('{}/{}/profiles_avg'.format(path, suffix))

        with Sync() as sync:
            if sync.comm.rank == 0:
                if not output_mid_path_yz.exists():
                    output_mid_path_yz.mkdir()
                if not output_mid_path_zx.exists():
                    output_mid_path_zx.mkdir()
                if not output_mid_path_xy.exists():
                    output_mid_path_xy.mkdir()
                # if not output_avg_path_yz.exists():
                #     output_avg_path_yz.mkdir()
                # if not output_avg_path_zx.exists():
                #     output_avg_path_zx.mkdir()
                # if not output_avg_path_xy.exists():
                #     output_avg_path_xy.mkdir()
                # if not output_path_avg.exists():
                #     output_path_avg.mkdir()

        post.visit_writes(slicepoints, plot_all)
        # post.visit_writes(slicepoints, yzmean, output=output_path_avg)