from configparser import ConfigParser
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

    filename = glob.glob("{}/options.cfg".format(sys.argv[1]))[0]
    print(filename)
    config = ConfigParser()
    config.read(str(filename))
    Rm = config.getfloat('parameters', 'Rm')

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
    if times_sorted[-1] > 5000:
        times_sorted = np.array(times_sorted)
        data_sorted = np.array(data_sorted)
        mask = (times_sorted >= 1000) & (times_sorted <= 10000)
        times_sorted = times_sorted[mask].tolist()
        data_sorted = data_sorted[mask].tolist()

    N = len(times_sorted)
    times_uniform = np.linspace(min(times_sorted), max(times_sorted), N)
    f = interpolate.interp1d(times_sorted, data_sorted)
    data_uniform = f(times_uniform)
    if 'floquet' in dir:
        idx = np.argmax(times_uniform >= 23.33)
        times_uniform = times_uniform[:idx + 1]
        data_uniform = data_uniform[:idx + 1]

    N = len(times_uniform)
    print(N)
    N = len(data_uniform)
    print(N)
    # sys.exit()

    if not 'floquet' in dir:
        window = signal.windows.hamming(len(data_uniform))
        wdata_uniform = window*data_uniform
    else:
        wdata_uniform = data_uniform
    # yf = fft(wdata_uniform)
    xf = fftfreq(N, d=times_uniform[1] - times_uniform[0])[:N//2]

    adjusted_0time = 23.33
    y_osc0 = np.sin(2*np.pi*times_uniform / adjusted_0time)
    yf_osc = fft(y_osc0)
    max_ind = np.argmax(2.0/N * np.abs(yf_osc[0:N//2]))
    plt.axvline(2/23.33, linestyle='dashed', color='black')
    # plt.axvline(1/23.33, linestyle='dotted', color='yellow')
    
    if 'floquet' in dir:
        wdata_uniform = data_uniform
        yf = fft(wdata_uniform)
        maxima_idx = np.where((np.diff(np.sign(np.diff(2.0/N * np.abs(yf[0:N//2])))) < 0))[0] + 1
        # plt.axvline(0.085419023875809, linestyle='dashed', color='orangered')
        # plt.axvline([xf[maxima_idx[0]]], linestyle='dashed', color='deeppink')
        # plt.axvline([xf[maxima_idx[1]]], linestyle='dashed', color='lime')
        # plt.axvline([xf[maxima_idx[2]]], linestyle='dashed', color='cyan')
        # plt.axvline([xf[maxima_idx[3]]], linestyle='dashed', color='darkviolet')

        T = 23.33
        plt.axvline(3/T, linestyle='dashed', color='deeppink')

        plt.axvline(5/T, linestyle='dashed', color='lime')

        # y_osc0 = np.sin(8*np.pi*times_uniform / adjusted_0time)
        # yf = fft(y_osc0)
        # max_ind = np.argmax(2.0/N * np.abs(yf[0:N//2]))
        plt.axvline(7/T, linestyle='dashed', color='cyan')

        plt.axvline(9/T, linestyle='dashed', color='darkviolet')

        cutoff = 12
        plt.scatter(xf[1:cutoff], 2.0/N * np.abs(yf[1:cutoff]), color='grey', alpha=1.0, zorder=10)

        print(maxima_idx)
        print(xf[1])
        print(xf[2])
        print(xf[maxima_idx])
        plt.xscale('linear')
        if '1e3' in dir:
            plt.ylabel('linear')
        plt.xticks([1/T, 3/T, 5/T, 7/T, 9/T], [r"$T^{-1}$", r"$3T^{-1}$", r"$5T^{-1}$", r"$7T^{-1}$", r"$9T^{-1}$"])
        # plt.xlim(0, 12/T)
        # plt.xticks([1/T, 3/T, 5/T, 7/T, 9/T], [r"$\frac{1}{T}$", r"$\frac{3}{T}$", r"$\frac{5}{T}$", r"$\frac{7}{T}$", r"$\frac{9}{T}$"])

    
    else:
        yf = fft(wdata_uniform)
        y_osc0 = np.sin(2*np.pi*times_uniform / 23.33 * 2)
        # yf = fft(y_osc0)
        # max_ind = np.argmax(2.0/N * np.abs(yf[0:N//2]))

        # plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]), color='black')
        # print(xf[np.argmax(2.0/N * np.abs(yf[0:N//2]))])
        # plt.axvline(0.08569953722085479, linestyle='dashed', color='orangered')
        # plt.axvline(0.21354756 / 2, linestyle='dashed', color='green')


        # plt.axvline(0.12812854, linestyle='dashed', color='deeppink')
        plt.axvline(0.12812854, linestyle='dashed', color='deeppink')

        plt.axvline(0.21354756, linestyle='dashed', color='lime')

        # y_osc0 = np.sin(8*np.pi*times_uniform / adjusted_0time)
        # yf = fft(y_osc0)
        # max_ind = np.argmax(2.0/N * np.abs(yf[0:N//2]))
        plt.axvline(0.29896658, linestyle='dashed', color='cyan')

        plt.axvline(0.38438561, linestyle='dashed', color='darkviolet')
        plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]), color='grey', alpha=0.8)
        plt.xscale('log')
        if '1e3' in dir:
            plt.ylabel('nonlinear')
        plt.xlim(0, 1)

    plt.yscale('log')
    # plt.xlabel('frequency')
    # plt.ylabel('')
    if 'floquet' in dir:
        plt.title(r"Rm = "+str(int(Rm)))
    # plt.ylabel(r"$E_y(f)$")
    plt.xlabel(r"$f$")


    plt.tight_layout()
    # plt.ylim(0.065, 0.1)
    figname = '{}oscillation_spectrum.png'.format(dir)
    plt.savefig(figname)
    print(figname)

    figname = '{}oscillation_spectrum.pdf'.format(dir)
    plt.savefig(figname)
    print(figname)
    

            # be_z = f['tasks']['be_z'][()].ravel()
            # be_x = f['tasks']['be_x'][()].ravel()
            # ke_y = f['tasks']['ke_y'][()].ravel()
            # ke_z = f['tasks']['ke_z'][()].ravel()
            # ke_x = f['tasks']['ke_x'][()].ravel()

# from configparser import ConfigParser
# try:
#     filename = glob.glob("{}/*.cfg".format(sys.argv[1]))[0]
#     config = ConfigParser()
#     config.read(str(filename))
# except:
#     filename = glob.glob("{}/{}/*.cfg".format(path, 'dnsNy4_T1e4'))[0]
#     config = ConfigParser()
#     config.read(str(filename))
#     print("failed to read config file. Please supply run suffix")
    # raise

# parent = config.get('parameters', 'suffix')
parent = sys.argv[1]
size = CW.size
rank = CW.rank 
targets = []
seeds = []

plot_spectrum(parent + '/')
