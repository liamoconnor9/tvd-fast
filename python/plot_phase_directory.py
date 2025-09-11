import os
import h5py
path = os.path.dirname(os.path.abspath(__file__))
import glob
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
data = []
import sys

dirname = "blend0"
suffices = glob.glob("{}/{}/*".format(path, dirname))
print(suffices)
golden_mean = (np.sqrt(5)-1.0)/2.0
ratio = 1.0
# ratio = 1.075
fig_width = 5 # in column width
# fig_width = 7.1 # in page width
fig_height = golden_mean * fig_width
fig_size =  [fig_width, fig_height]
# plt.figure(figsize=fig_size)
plt.rcParams.update({'font.size': 16})
# fig = plt.figure(figsize=(6, 6))


def plot_23_phase(dir, prefix='ke'):

    files = glob.glob("{}/scalars/*.h5".format(dir))
    last_index = len(files)
    first = True
    for file in files:
        if "_s{}.h5".format(last_index) in file:
            print("success")
            # continue
        labels = {}
        with h5py.File(file, "r") as f:
            Nmodes = 4
            ke_modes = []
            for ki in range(Nmodes):
                ke_modes.append(f['tasks']['{}_mode{}'.format(prefix, ki + 1)][()].ravel())
                labels['{}_mode{}'.format(prefix, ki + 1)] = ke_modes[-1]

        plt.plot(labels['ke_mode2'], labels['ke_mode3'], color='grey', alpha=0.3)

for ind, suffix in enumerate(suffices):
    plot_23_phase(suffix + '/')

plt.xlabel("mode 2")
plt.ylabel("mode 3")
plt.gca().set_aspect('equal')
plt.xlim(0, 0.07)
plt.ylim(0, 0.07)

plt.legend()
plt.tight_layout()
figname = path + '/{}/phase_23.png'.format(dirname)
plt.savefig(figname)
print(figname)
