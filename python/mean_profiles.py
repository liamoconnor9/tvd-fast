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


def yzmean(filename, start, count):
    """Save plot of specified tasks for given range of analysis writes."""
    if "_s{}.h5".format(last_index) in filename:
        print('success')
        return
    
    minimum_index = 20
    for index in range(minimum_index):
        if "_s{}.h5".format(index) in filename:
            print('bypassing file {}'.format(filename))
            return


    global sum_sum_vec, x, tasks
    # Plot settings
    normal_dir = 'z'
    scale = 2.5
    dpi = 100
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: 'mid_{:06}.png'.format(write)
    # Layout
    # if (round(ary) > 1):
    nrows, ncols = 2, 3


    # Plot writes
    with h5py.File(filename, mode='r') as file:
        for key in file['scales'].keys():
            if 'x_hash' in key:
                x = file['scales'][key][()]
        sum_vec = [np.zeros_like(file['tasks'][tasks[0]][()][0, 0, 0, :]) for task in tasks]
        sum_vec.append(0)
        for index in range(start, start+count):
            sim_time = file['scales/sim_time'][index]
            timestep = file['scales/timestep'][index]
            sum_vec[-1] += 1
            for n, task in enumerate(tasks):
                sum_vec[n] += file['tasks'][task][()][index, 0, 0, :]

    if len(sum_sum_vec) == 0:
        sum_sum_vec = sum_vec
    else:
        for i in range(len(sum_vec)):
            sum_sum_vec[i] += sum_vec[i]
        


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
    filename = "{}/{}/mri_options.cfg".format(path, suffix)
    config = ConfigParser()
    config.read(str(filename))

    global sum_sum_vec, x, ar, ary, arz, last_index, tasks
    tasks = ['vy_avg', 'by_avg', 'jy_avg', 'vz_avg', 'bz_avg', 'jz_avg']
    sum_sum_vec = []
    Ly = eval(config.get('parameters','Ly'))
    Lz = eval(config.get('parameters','Lz'))
    Lx = eval(config.get('parameters','Lx'))

    ary = Ly / Lx
    arz = Lz / Lx 

    data_dirs = glob("{}/{}/data/*/".format(path, suffix))
    for data_dir in data_dirs:
        index = int(data_dir[:-1].split("/")[-1])
        slicepoints = glob("{}/{}/data/{}/slicepoints/*.h5".format(path, suffix, index))
        last_index = len(slicepoints)

        post.visit_writes(slicepoints, yzmean)

    write_dict = {}    
    total_time = sum_sum_vec[-1]
    for n, element in enumerate(sum_sum_vec):
        element /= total_time
        if n < len(tasks):
            write_dict[tasks[n]] = element

    write_dict['x'] = x
    np.save("{}/{}/data/1/write_dict.npy".format(path, suffix), write_dict)