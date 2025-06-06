from scipy import stats
from scipy.optimize import root, root_scalar, minimize_scalar
import numpy as np
import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
import h5py
import dedalus.public as d3
from mpi4py import MPI
CW = MPI.COMM_WORLD
rank = CW.Get_rank()
sys.path.append("..") # Adds higher directory to python modules path.
from vp_bvp_func import *
from glob import glob
import logging
logger = logging.getLogger(__name__)
from docopt import docopt
from pathlib import Path
from configparser import ConfigParser
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from run_sim import run_sim

filename = path + "/options.cfg"
ky = float(sys.argv[1])
# Rm = 36.86
global config
config = ConfigParser()
config.read(str(filename))
Rm = config.getfloat('parameters', 'Rm')
# logger.info(config.items('parameters'))

def write_output(output, mode):
    with open('root_output0.txt',mode) as file:
        file.write(output)
        file.write('\n')

def growth_rate(ky, Rm):
    global config
    logger.info("RUNNING KINEMATIC EXP. FOR Rm, ky = {}, {}".format(Rm, ky))
    time_lst, be_y_lst = run_sim(ky, Rm, config)
    times_data = np.array(time_lst)
    data_data = np.array(be_y_lst)
    if rank == 0:
        x, y = times_data, np.log(data_data)
        plt.plot(x, y)
        plt.title('ky={}'.format(ky))
        plt.xlabel('time')
        plt.ylabel('log(azimuthal BE)')
        figname = 'ky{}'.format(str(ky).replace('.', 'p'))
        plt.savefig("{}/traces/{}.png".format(path, figname))
        plt.close()

        slope, intercept, r, p, std_err = stats.linregress(x, y)
        msg = "{}, {}, {}".format(ky, slope, std_err)
        write_output(msg, 'a')
    else:
        slope = np.inf
    CW.Barrier()
    slope = CW.bcast(slope, root=0)
    CW.Barrier()
    return 

growth_rate(ky, Rm)
# write_output('ky, growth_rate, std_err', 'w')
