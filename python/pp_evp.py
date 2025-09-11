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

# ary = Ly / Lx
# arz = Lz / Lx 

# data_dirs = glob("{}/{}/data/*/".format(path, suffix))
# for data_dir in data_dirs:
#     index = int(data_dir[:-1].split("/")[-1])
#     slicepoints = glob("{}/{}/data/{}/slicepoints/*.h5".format(path, suffix, index))
#     last_index = len(slicepoints)

#     post.visit_writes(slicepoints, yzmean)

# write_dict = {}    
# total_time = sum_sum_vec[-1]
# for n, element in enumerate(sum_sum_vec):
#     element /= total_time
#     if n < len(tasks):
#         write_dict[tasks[n]] = element

data = np.load("{}/{}/data/1/write_dict.npy".format(path, suffix), allow_pickle=True)[()]
x = data['x']

fig, axes = plt.subplots(3, 3, sharex=True, sharey=False, layout='constrained')
for n, task in enumerate(tasks):
    # Build subfigure axes
    i, j = divmod(n, 3)
    dset_vec = data[task]
    axes[i, j].plot(x, dset_vec, linewidth=3, color='purple')
    axes[i, j].set_title(task)
    if i == 2:
        axes[i, j].set_xlabel('x')
    axes[i, j].set_xlim(-0.5, 0.5)

# # Add time title
# # Save figure
figname = "{}/{}/data/1/mean_profiles.png".format(path, suffix)
plt.savefig(figname, dpi=100)
print(figname)

