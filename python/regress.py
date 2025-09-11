from docopt import docopt
from configparser import ConfigParser
from pathlib import Path
import numpy as np
import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
import h5py
import dedalus.public as d3
from mpi4py import MPI
CW = MPI.COMM_WORLD
import logging
import pathlib
from glob import glob
logger = logging.getLogger(__name__)
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

try:    
    args = docopt(__doc__)
    filename = Path(args['<config_file>'])
    seed = 1
    # seed = int(args['<seed>'])
except:
    logger.warning('invalid config file supplied. using default... w/ seed = 1')
    filename = path + "/mri_options.cfg"
    seed = 1

config = ConfigParser()
config.read(str(filename))

logger.info(config.items('parameters'))

scale = config.getfloat('parameters','scale')
logger_cadence = config.getint('parameters','logger_cadence')

suffix = eval(config.get('parameters', 'suffix'))




def get_periods(dir, task, logscale=False):
    files = glob("{}scalars/*.h5".format(dir))
    files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    last_index = len(files)
    times_data = []
    data_data = []
    for file in files:
        # print(file)
        if "_s{}.h5".format(last_index) in file:
            print("success")
            continue
        with h5py.File(file, "r") as f:
            times_data += f['scales']['sim_time'][()].tolist()
            data_data += f['tasks'][task][()].ravel().tolist()

    plt.plot(times_data, data_data, color='purple')
    peaks, _ = find_peaks(data_data)
    periods = [times_data[peaks[indy + 1]] - times_data[peaks[indy]] for indy in range(len(peaks) - 1)]
    return periods

def get_vals(dir, task, time_pt, period):
    files = glob("{}scalars/*.h5".format(dir))
    files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    last_index = len(files)
    times_data = []
    data_data = []
    for file in files:
        # print(file)
        if "_s{}.h5".format(last_index) in file:
            print("success")
            continue
        with h5py.File(file, "r") as f:
            times_data += f['scales']['sim_time'][()].tolist()
            data_data += f['tasks'][task][()].ravel().tolist()
    
    val0 = np.interp(time_pt, times_data, data_data)
    val1 = np.interp(time_pt+period, times_data, data_data)

    return (val0, val1)

from configparser import ConfigParser
try:
    filename = glob("{}/*.cfg".format(sys.argv[1]))[0]
    config = ConfigParser()
    config.read(str(filename))
    Pm = config.getfloat("parameters", "Pm")
    # try:
    #     Ro = config.getfloat("parameters", "Ro")
    # except:
    #     t_vec = eval(eval(config.get("parameters", "t_vec")))
    #     Ro_vec = (eval(config.get("parameters", "Ro_vec")))
    #     Ro = (Ro_vec[0], Ro_vec[-1])
    # icContinuation = config.getfloat("parameters", "Ro")

    nu = config.getfloat("parameters", "nu")
except:
    print("failed to read config file. Please supply run suffix")
    raise

parent = sys.argv[1] + '/data/'
size = CW.size
rank = CW.rank 
targets = []
seeds = []
for folder in glob(parent + "*/"):
    folder = folder.replace(parent, "")
    folder = folder.replace("/", "")
    if folder.isdigit():
        targets.append(folder)
        seeds.append(int(folder))

N = len(list(range(rank, len(seeds), size)))
for i in range(rank, len(seeds), size):


    periods = get_periods(parent + targets[i] + '/', 'udiff')
    pruneBeginning = True
    # print('standard deviations')
    prune_indy = 0
    for j in range(1,int(len(periods)/2)):
        if np.std(periods[j:]) < 1e-3:
            prune_indy = j
            break
    period = np.mean(periods[prune_indy:])
    def gr_from_vals(args):
        return args[1] / args[0]
    print(gr_from_vals(get_vals(parent + targets[i] + '/', 'be_y', 75, period)))
    print(gr_from_vals(get_vals(parent + targets[i] + '/', 'be_y', 76, period)))
    print(gr_from_vals(get_vals(parent + targets[i] + '/', 'be_y', 77, period)))

    # print(period)
    sys.exit()