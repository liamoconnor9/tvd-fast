import numpy as np
# import dedalus.public as d3
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from mpi4py import MPI
CW = MPI.COMM_WORLD

from configparser import ConfigParser
try:
    filename = glob.glob("{}/*.cfg".format(sys.argv[1]))[0]
    config = ConfigParser()
    config.read(str(filename))
except:
    print("failed to read config file. Please supply run suffix")
    raise

parent = sys.argv[1]
size = CW.size
rank = CW.rank 
# targets = []
# seeds = []
# for folder in glob.glob(parent + "*/"):
#     folder = folder.replace(parent, "")
#     folder = folder.replace("/", "")
#     if folder.isdigit():
#         targets.append(folder)
#         seeds.append(int(folder))

# N = len(list(range(rank, len(seeds), size)))


dir = parent + '/'
# task = 'udiff'
task = 'projA'
logscale = True

files = glob.glob("{}scalars/*.h5".format(dir))
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



from scipy.signal import find_peaks

if 'proj' in task:
    data_data = 1 - np.array(data_data)
    peaks_full, _ = find_peaks(-np.array(data_data), height=-0.1)
else:
    peaks_full, _ = find_peaks(-np.array(data_data))
# print(peaks_full)
times_data = np.array(times_data)
t_f = times_data[peaks_full]
# print(np.diff(t_f))
periods = np.diff(t_f)
print('t_f = {}'.format(t_f))
print('mean period = {}'.format(np.mean(periods)))
print('periods = {}'.format(periods))

plt.hist(periods)
plt.xlabel("periods")
plt.ylabel("probability")
# if logscale:
#     plt.yscale('log')
# plt.yscale('log')
figname = "{}{}.png".format(dir, 'periods_hist')
plt.savefig(figname)
print(figname)
plt.close()

# for i in range(rank, len(seeds), size):
#     try:
#         plot_task(parent + targets[i] + '/', 'udiff', logscale=True)
#     except:
#         print('failed to plot udiff')