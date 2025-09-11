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

# Apply a window function
def plot_spectrum(dir):
    plt.figure(figsize=(4, 3))
    files = glob.glob("{}scalars/*.h5".format(dir))
    # last_index = len(files)
    times_lst = []
    data_lst = []
    for file in files:
        # if "_s{}.h5".format(last_index) in file:
        #     print("success")
        #     continue

        with h5py.File(file, "r") as f:

            sim_times = f['scales']['sim_time'][()]
            # be_y = f['tasks']['be_y'][()].ravel()
            be_y = f['tasks']['be_y'][()].ravel()

            times_lst += sim_times.tolist()
            data_lst += be_y.tolist()
            # print(times_lst)
            # sys.exit()

    times_sorted = [x for x, _ in sorted(zip(times_lst, data_lst))]
    data_sorted = [x for _, x in sorted(zip(times_lst, data_lst))]
    # print(data_sorted)
    N = len(times_sorted)
    times_uniform = np.linspace(min(times_sorted), max(times_sorted), N)
    f = interpolate.interp1d(times_sorted, data_sorted)
    data_uniform = f(times_uniform)
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
        # plt.axvline(0.085419023875809, linestyle='dashed', color='orangered')
        plt.axvline([xf[maxima_idx[0]]], linestyle='dashed', color='deeppink')
        plt.axvline([xf[maxima_idx[1]]], linestyle='dashed', color='lime')
        plt.axvline([xf[maxima_idx[2]]], linestyle='dashed', color='cyan')
        plt.axvline([xf[maxima_idx[3]]], linestyle='dashed', color='darkviolet')
        print(maxima_idx)
        print(xf[1])
        print(xf[2])
        print(xf[maxima_idx])

    
    truncate = -1
    if not 'floquet' in dir:
        yf = fft(wdata_uniform)
        plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]), color='grey')
        y_osc0 = np.sin(2*np.pi*times_uniform / 150)
        yf = fft(y_osc0)
        max_ind = np.argmax(2.0/N * np.abs(yf[0:N//2]))

        # plt.axvline(0.085419023875809, linestyle='dashed', color='orangered')


        plt.axvline(0.12812854, linestyle='dashed', color='deeppink')

        plt.axvline(0.21354756, linestyle='dashed', color='lime')

        # y_osc0 = np.sin(8*np.pi*times_uniform / adjusted_0time)
        # yf = fft(y_osc0)
        # max_ind = np.argmax(2.0/N * np.abs(yf[0:N//2]))
        plt.axvline(0.29896658, linestyle='dashed', color='cyan')

        plt.axvline(0.38438561, linestyle='dashed', color='darkviolet')

    adjusted_0time = 23.33
    y_osc0 = np.sin(4*np.pi*times_uniform / adjusted_0time)
    yf = fft(y_osc0)
    max_ind = np.argmax(2.0/N * np.abs(yf[0:N//2]))
    plt.axvline(xf[max_ind], linestyle='dashed', color='orangered')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(0, 1)
    plt.tight_layout()
    # plt.ylim(0.065, 0.1)
    figname = '{}oscillation_spectrum.png'.format(dir)
    plt.savefig(figname)
    print(figname)
    

            # be_z = f['tasks']['be_z'][()].ravel()
            # be_x = f['tasks']['be_x'][()].ravel()
            # ke_y = f['tasks']['ke_y'][()].ravel()
            # ke_z = f['tasks']['ke_z'][()].ravel()
            # ke_x = f['tasks']['ke_x'][()].ravel()

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
parent = sys.argv[1]
size = CW.size
rank = CW.rank 
targets = []
seeds = []

plot_spectrum(parent + '/')
