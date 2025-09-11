import numpy as np
# import dedalus.public as d3
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from mpi4py import MPI
CW = MPI.COMM_WORLD
from scipy import interpolate
from scipy.fft import fft, fftfreq
import scipy.signal as signal
from scipy.signal import savgol_filter
import os
path = os.path.dirname(os.path.abspath(__file__))

plt.figure(figsize=(8, 6))
curve_ind = 0
colors = ['purple', 'navy']
linestyles = ['solid', 'solid']

# Apply a window function
def plot_spectrum(dir, label):
    global curve_ind, colors, linestyles
    files = glob.glob("{}scalars/*.h5".format(dir))
    # last_index = len(files)
    times_lst = []
    data_lst = []
    for file in files:
        with h5py.File(file, "r") as f:
            sim_times = f['scales']['sim_time'][()]
            be_y = f['tasks']['be_y'][()].ravel()
            times_lst += sim_times.tolist()
            data_lst += be_y.tolist()

    times_sorted = [x for x, _ in sorted(zip(times_lst, data_lst))]
    data_sorted = [x for _, x in sorted(zip(times_lst, data_lst))]
    # print(data_sorted)
    N = len(times_sorted)
    times_uniform = np.linspace(min(times_sorted), max(times_sorted), N)
    f = interpolate.interp1d(times_sorted, data_sorted)
    data_uniform = f(times_uniform)
    # data_uniform = np.sin(times_uniform * 2 * np.pi / 100)
    if 'floquet' in dir:
        idx = np.argmax(times_uniform >= 23.33)
        times_uniform = times_uniform[:idx + 1]
        data_uniform = data_uniform[:idx + 1]

    # Number of sample points
    window = signal.windows.hamming(len(data_uniform))
    wdata_uniform = window*data_uniform
    # yf = fft(wdata_uniform)
    xf = fftfreq(N, d=times_uniform[1] - times_uniform[0])[:N//2]
    if 'floquet' in dir:
        wdata_uniform = data_uniform
        yf = fft(wdata_uniform)
        maxima_idx = np.where((np.diff(np.sign(np.diff(2.0/N * np.abs(yf[0:N//2])))) < 0))[0] + 1
        plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]), color='grey')
    
    truncate = -1
    if not 'floquet' in dir:
        yf = fft(wdata_uniform)
        plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]), linewidth=0.5, linestyle=linestyles[curve_ind], color=colors[curve_ind], label=label)
    curve_ind += 1

from configparser import ConfigParser
try:
    filename = glob.glob("{}/*.cfg".format(sys.argv[1]))[0]
    config = ConfigParser()
    config.read(str(filename))
except:
    filename = glob.glob("{}/{}/*.cfg".format(path, 'dnsNy4_T1e4'))[0]
    config = ConfigParser()
    config.read(str(filename))
    print("failed to read config file. Please supply run suffix")
    # raise

# parent = config.get('parameters', 'suffix')
for parent in sys.argv[1:]:
    print(parent)
    filename = glob.glob("{}/*.cfg".format(parent))[0]
    config = ConfigParser()
    config.read(str(filename))
    Rm = int(config.getfloat('parameters', 'Rm'))
    label = "Rm = {}".format(Rm)
    plot_spectrum(parent + '/', label)
parent = sys.argv[1]

plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlim(2e-4, 1)
plt.tight_layout()
# plt.ylim(0.065, 0.1)
figname = '{}/combined_spectrum.png'.format(parent)
plt.savefig(figname)
print(figname)
