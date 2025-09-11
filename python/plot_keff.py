import matplotlib
matplotlib.use('agg')
from scipy import stats
import h5py
import glob
import sys
import matplotlib.pyplot as plt
import numpy as np
import os
path = os.path.dirname(os.path.abspath(__file__))
from configparser import ConfigParser
def analyze_suffix(parent):
    print(parent)
    try:
        filename = glob.glob("{}/{}/*.cfg".format(path, parent))[0]
        config = ConfigParser()
        config.read(str(filename))
        Re = config.getfloat("parameters", "Re")
        Rm = config.getfloat("parameters", "Rm")
    except:
        print("failed to read config file. Please supply run suffix: " + parent)
        raise
    dir = parent + '/'
    task = 'keff'
    files = glob.glob("{}/{}scalars/*.h5".format(path, dir))
    files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    # print('here')
    times_data = []
    data_data = []
    for file in files:
        with h5py.File(file, "r") as f:
            times_data += f['scales']['sim_time'][()].tolist()
            data_data += f['tasks'][task][()].ravel().tolist()

    data_data = np.array(data_data)
    times_data = np.array(times_data)
    if times_data[-1] < 23.33:
        print("config not run long enough")
        return None
    # Ncutoff = len([time for time in times_data if time <= times_data[-1] - 23.33])
    # times_data = times_data[Ncutoff:]
    # data_data = data_data[Ncutoff:]
    return (Rm, np.mean(data_data))

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(6, 4))
markersize1 = 100
markersize2 = 100
marker1 = 'o'
marker2 = '^'
alpha = 0.5
suffices = []

with open('{}/suffices_keff_Ly22.txt'.format(path), 'r') as file:
    for i in range(0, 90):
        temp = file.readline().replace('\n', '')
        if len(temp) == 0:
            break
        suffices.append(temp)

rm = []
gr = []
for i, suffix in enumerate(suffices):
    if i == 0:
        continue
    try:
        val = analyze_suffix(suffix)
    except Exception as e:
        print('failed for suffix = {}'.format(suffix))
        print(e)
        continue
    
    if val == None:
        continue
    else:
        rm.append(val[0])
        gr.append(val[1])

gr = np.array(gr)
rm = np.array(rm)
keff = gr**2 / rm
plt.scatter(rm, keff, color='magenta', label=r"$k_y=0.29$", marker=marker1, s=markersize1, alpha=alpha)

suffices = []
with open('{}/suffices_keff_Ly2pi.txt'.format(path), 'r') as file:
    for i in range(0, 90):
        temp = file.readline().replace('\n', '')
        if len(temp) == 0:
            break
        suffices.append(temp)

rm = []
gr = []
for suffix in suffices:
    if i == 0:
        continue
    try:
        val = analyze_suffix(suffix)
    except Exception as e:
        print('failed for suffix = {}'.format(suffix))
        print(e)
        continue

    if val == None:
        continue
    else:
        rm.append(val[0])
        gr.append(val[1])

gr = np.array(gr)
rm = np.array(rm)
keff = gr**2 / rm
plt.scatter(rm, keff, color='cyan', label=r"$k_y=1.0$", marker=marker2, s=markersize2, alpha=alpha)

plt.legend(frameon=False)
plt.xscale('log')
plt.xlabel(r'$\rm{Rm}$')
plt.ylabel(r'$k_{\rm{eff}}^2 \, / \, \rm{Rm}$')
plt.tight_layout()

figname = "{}/keff.png".format(path)
plt.savefig(figname)
print(figname)

figname = "{}/keff.pdf".format(path)
plt.savefig(figname)
print(figname)