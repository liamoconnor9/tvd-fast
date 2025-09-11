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
from mpi4py import MPI
CW = MPI.COMM_WORLD

def plot_plane(filename, start, count, output, normal_dir, tag):
    """Save plot of specified tasks for given range of analysis writes."""

    if normal_dir == 'x':
        image = plot_tools.Box(2, 2 * ary / arz)
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
    dpi = 200
    title_func = lambda sim_time: 't = {:.1f}'.format(sim_time)
    savename_func = lambda write: '{}_{:06}.png'.format(tag, write)
    nrows, ncols = 1, 3
    tasks = ['vy_midy', 'b.grad_b{}_midy'.format(tag[-1]), 'b.grad_b{}_avgy'.format(tag[-1])]

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
            data_slices = (index, ) + data_slices_tail
            for n, task in enumerate(tasks):
                # Build subfigure axes
                i, j = divmod(n, ncols)
                axes = mfig.add_axes(i, j, [0, 0, 1, 1])
                dset = file['tasks'][task]
                if 'b' in task:
                    pl = plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes, title=task, even_scale=True, cmap="PiYG")
                    cb_ax = pl[1]
                    # cb_ax.set_xticks([vmin_b, 0, vmax_b])
                else:
                    pl = plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes, title=r"$u_y$", even_scale=True, clim=(vmin_u, vmax_u))
                    cb_ax = pl[1]
                    cb_ax.set_xticks([vmin_u, 0, vmax_u])
            # print()
            # sys.exit()
            # im = plt.gca().images        
            # cb = im[-1].colorbar   
            # print(cb)

            # Add time title
            title = title_func(file['scales/sim_time'][index])
            title_height = 1 - 0.5 * mfig.margin.top / mfig.fig.y
            fig.suptitle(title, x=0.42, y=title_height, ha='left')
            # Save figure
            savename = savename_func(file['scales/write_number'][index])
            savepath = output.joinpath(savename)
            if (index % 1 == 0):
                fig.savefig(str(savepath), dpi=dpi)
            fig.clear()
    plt.close(fig)

def plot_all(filename, start, count):
    if "_s{}.h5".format(last_index) in filename:
        print('success')
    plot_plane(filename, start, count, output_mid_path_bdotgradby, 'y', 'bdotgradby')
    plot_plane(filename, start, count, output_mid_path_bdotgradbz, 'y', 'bdotgradbz')
    plot_plane(filename, start, count, output_mid_path_bdotgradbx, 'y', 'bdotgradbx')


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

    global ar, ary, arz, last_index, isHydro, solver_name
    Ly = eval(config.get('parameters','Ly'))
    Lz = eval(config.get('parameters','Lz'))
    Lx = eval(config.get('parameters','Lx'))
    isHydro = config.getboolean('parameters','isHydro')
    solver_name = eval(config.get('parameters','SOLVER'))

    ary = Ly / Lx
    arz = Lz / Lx 

    slicepoints = glob("{}/{}/slicepoints/*.h5".format(path, suffix))
    last_index = len(slicepoints)

    # Create output directory if needed
    output_mid_path_bdotgradby=pathlib.Path('{}/{}/bdotgradby'.format(path, suffix))
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_mid_path_bdotgradby.exists():
                output_mid_path_bdotgradby.mkdir()

    output_mid_path_bdotgradbz=pathlib.Path('{}/{}/bdotgradbz'.format(path, suffix))
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_mid_path_bdotgradbz.exists():
                output_mid_path_bdotgradbz.mkdir()
                
    output_mid_path_bdotgradbx=pathlib.Path('{}/{}/bdotgradbx'.format(path, suffix))
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_mid_path_bdotgradbx.exists():
                output_mid_path_bdotgradbx.mkdir()

    post.visit_writes(slicepoints, plot_all)