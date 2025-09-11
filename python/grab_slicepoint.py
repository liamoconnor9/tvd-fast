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

t1 = 5.0
t2 = 10.8325
# quarter period is 5.8325
vy_t1 = np.zeros((1, 256, 64))
vy_t2 = np.zeros((1, 256, 64))
vz_t1 = np.zeros((1, 256, 64))
vz_t2 = np.zeros((1, 256, 64))
vx_t1 = np.zeros((1, 256, 64))
vx_t2 = np.zeros((1, 256, 64))
by_t1 = np.zeros((1, 256, 64))
by_t2 = np.zeros((1, 256, 64))
done1 = done2 = False

def full_dot(mat1, mat2):
    return np.sum(np.multiply(mat1, mat2))

def plot_plane(filename, start, count, normal_dir, tag):
    """Save plot of specified tasks for given range of analysis writes."""
    global done1, done2, t1, t2, vy_t1, vy_t2, vz_t1, vz_t2, vx_t1, vx_t2, by_t1, by_t2

    if done1 and done2:
        return
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

    # Plot writes
    with h5py.File(filename, mode='r') as file:
        for index in range(start, start+count):
            data_slices = (index, ) + data_slices_tail
            if np.abs(file['scales']['sim_time'][index] - t1) < 0.001:
                vy_t1 = file['tasks']['vy_midy'][index][()]
                vz_t1 = file['tasks']['vz_midy'][index][()]
                vx_t1 = file['tasks']['vx_midy'][index][()]
                by_t1 = file['tasks']['by_midy'][index][()]
                done1 = True

            if np.abs(file['scales']['sim_time'][index] - t2) < 0.001:
                vy_t2 = file['tasks']['vy_midy'][index][()]
                vz_t2 = file['tasks']['vz_midy'][index][()]
                vx_t2 = file['tasks']['vx_midy'][index][()]
                by_t2 = file['tasks']['by_midy'][index][()]
                done2 = True

        if done1 and done2:
            print('done')
            print(by_t1[0, 45, 48])
            print(by_t2[0, 45, 48])
            print(t1)
            print(t2)
            return
def plot_all(filename, start, count):
    plot_plane(filename, start, count, 'y', 'mid')

if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync
    suffix = 'kin_2p6d_Rm1500_RSTRT2'

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
        post.visit_writes(slicepoints, plot_all)
        # post.visit_writes(slicepoints, yzmean, output=output_path_avg)
    
    np.savez('{}/snaps_saved'.format(path), vy_t1=vy_t1, vy_t2=vy_t2, by_t1=by_t1, by_t2=by_t2, t1=t1, t2=t2)
    # global done1, done2, t1, t2, vy_t1, vy_t2, by_t1, by_t2
