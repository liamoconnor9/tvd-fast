import numpy as np
# import dedalus.public as d3
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from mpi4py import MPI
CW = MPI.COMM_WORLD
import os
path = os.path.dirname(os.path.abspath(__file__))
from configparser import ConfigParser
from scipy import stats
from scipy.signal import find_peaks, peak_prominences
from scipy import interpolate

try:
    filename = glob.glob("{}/*.cfg".format(sys.argv[1]))[0]
    config = ConfigParser()
    config.read(str(filename))
    Re = config.getfloat("parameters", "Re")
    Rm = config.getfloat("parameters", "Rm")
except:
    print("failed to read config file. Please supply run suffix")
    raise

# kin_2p6d_Rm2500_RSTRT2
parent = sys.argv[1]
dir = parent + '/'
task = 'be_y'
logscale = True

files = glob.glob("{}scalars/*.h5".format(dir))
files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
last_index = len(files)
times_data = []
data_data = []
for index, file in enumerate(files):
    # print(file)
    # if "_s{}.h5".format(last_index) in file:
    #     print("success")
    #     continue
    with h5py.File(file, "r") as f:
        times_data += f['scales']['sim_time'][()].tolist()
        data_data += f['tasks'][task][()].ravel().tolist()

# find_peaks(data_data, height=None, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)

data_data = np.array(data_data)
times_data = np.array(times_data)

delta_data = [np.abs(data_data[i+1] - data_data[i]) for i in range(len(data_data) - 1)]
deltas_sorted = sorted(delta_data)[::-1]
# print(deltas_sorted[:15])
avg_delta_data = np.mean(delta_data)
last_peak = np.nan
# for i in range(len(delta_data)):
end_time = times_data[-1]
delta_min = 1
for i in range(len(delta_data) - 1, 1, -1):
    if delta_data[i] > delta_min:
        end_time = times_data[i - 2]
        break
print(end_time)
# sys.exit()

ln_data = np.log(data_data)
period = 23.33
first_time = end_time - 2*period
second_time = end_time - period
t_pts = [first_time, second_time, end_time]

f = interpolate.interp1d(times_data, ln_data)
y_pts = f(t_pts)

slope1 = (y_pts[1] - y_pts[0]) / (t_pts[1] - t_pts[0])
slope2 = (y_pts[2] - y_pts[1]) / (t_pts[2] - t_pts[1])
error = abs(slope1 - slope2) / min(abs(slope1), abs(slope2))
print("slope 1 = {}".format(slope1))
print("slope 2 = {}".format(slope2))
print("error = {}".format(error))
