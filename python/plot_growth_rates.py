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
    try:
        filename = glob.glob("{}/{}/*.cfg".format(path, parent))[0]
        config = ConfigParser()
        config.read(str(filename))
        Re = config.getfloat("parameters", "Re")
        Rm = config.getfloat("parameters", "Rm")
    except:
        print("failed to read config file. Please supply run suffix")
        raise
    dir = parent + '/'
    task = 'be_y'
    files = glob.glob("{}/{}scalars/*.h5".format(path, dir))
    files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    times_data = []
    data_data = []
    for file in files:
        with h5py.File(file, "r") as f:
            times_data += f['scales']['sim_time'][()].tolist()
            data_data += f['tasks'][task][()].ravel().tolist()

    data_data = np.array(data_data)
    times_data = np.array(times_data)
    # print(data_data)
    
    x, y = times_data, np.log(data_data)
    slope, intercept, r, p, std_err = stats.linregress(x, y)
    msg = "{}; rm={}; T={}; slope={}; rsqrd={}".format(parent, Rm, times_data[-1], slope, r**2)
    if r**2 > 0.99:
        print(msg)
        return (Rm, slope / 2)
    else:
        print(parent)
        return None


plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(6, 4))
markersize1 = 100
markersize2 = 100
marker1 = 'o'
marker2 = '^'
alpha = 0.5

suffices = []
with open('{}/critical_suffices.txt'.format(path), 'r') as file:
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
    except:
        print('failed for suffix = {}'.format(suffix))
        continue
    if val == None:
        continue
    else:
        rm.append(val[0])
        gr.append(val[1])
plt.scatter(rm, gr, color='magenta', label=r"$k_y=0.29$", marker=marker1, s=markersize1, alpha=alpha)

suffices = []
with open('{}/suffices_rm.txt'.format(path), 'r') as file:
    for i in range(0, 90):
        temp = file.readline().replace('\n', '')
        if len(temp) == 0:
            break
        suffices.append(temp)

rm = []
gr = []
for suffix in suffices:
    try:
        val = analyze_suffix(suffix)
    except:
        print('failed for suffix = {}'.format(suffix))
    if val == None:
        continue
    else:
        rm.append(val[0])
        gr.append(val[1])
plt.scatter(rm, gr, color='cyan', label=r"$k_y=1.0$", marker=marker2, s=markersize2, alpha=alpha)

plt.legend(frameon=False)
plt.xscale('log')
plt.xlabel(r'$\rm{Rm}$')
plt.ylabel('growth rate')
plt.tight_layout()

figname = "{}/grrm_auto3.png".format(path)
plt.savefig(figname)
print(figname)

figname = "{}/grrm_auto3.pdf".format(path)
plt.savefig(figname)
print(figname)