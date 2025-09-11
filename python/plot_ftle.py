from scipy import stats
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
import os
path = os.path.dirname(os.path.abspath(__file__))
from matplotlib.gridspec import GridSpec
import h5py

ftle1 = np.zeros((512, 128))
ftle2 = np.zeros((512, 128))

if True:
    loaded_NW_1 = np.load('/home/x-loconnor/mhd/d2/quadrants/lya_NW_1_0001.npz')
    N = np.shape(loaded_NW_1['lya'])[0]
    dt = 1e-2
    times = np.array([i*dt for i in range(N)])
    reshaped_NW_1 = np.reshape(loaded_NW_1['lya'], (N, 256, 64))

    loaded_SW_1 = np.load('/home/x-loconnor/mhd/d2/quadrants/lya_SW_1_0001.npz')
    reshaped_SW_1 = np.reshape(loaded_SW_1['lya'], (N, 256, 64))

    loaded_SE_1 = np.load('/home/x-loconnor/mhd/d2/quadrants/lya_SE_1_0001.npz')
    reshaped_SE_1 = np.reshape(loaded_SE_1['lya'], (N, 256, 64))

    loaded_NE_1 = np.load('/home/x-loconnor/mhd/d2/quadrants/lya_NE_1_0001.npz')
    reshaped_NE_1 = np.reshape(loaded_NE_1['lya'], (N, 256, 64))

    reshaped = reshaped_SW_1
    offset_row, offset_col = 0, 0
    for row in range(256):
        for col in range(64):
            data_vec = reshaped[:, row, col]
            x, y = times, data_vec
            slope, intercept, r, p, std_err = stats.linregress(x, y)
            ftle1[row + offset_row, col + offset_col] = slope

    reshaped = reshaped_SE_1
    offset_row, offset_col = 0, 64
    for row in range(256):
        for col in range(64):
            data_vec = reshaped[:, row, col]
            x, y = times, data_vec
            slope, intercept, r, p, std_err = stats.linregress(x, y)
            ftle1[row + offset_row, col + offset_col] = slope

    reshaped = reshaped_NW_1
    offset_row, offset_col = 256, 0
    for row in range(256):
        for col in range(64):
            data_vec = reshaped[:, row, col]
            x, y = times, data_vec
            slope, intercept, r, p, std_err = stats.linregress(x, y)
            ftle1[row + offset_row, col + offset_col] = slope

    reshaped = reshaped_NE_1
    offset_row, offset_col = 256, 64
    for row in range(256):
        for col in range(64):
            data_vec = reshaped[:, row, col]
            x, y = times, data_vec
            slope, intercept, r, p, std_err = stats.linregress(x, y)
            ftle1[row + offset_row, col + offset_col] = slope

    loaded_NW_33 = np.load('/home/x-loconnor/mhd/d2/quadrants/lya_NW_33_0033.npz')
    N = np.shape(loaded_NW_33['lya'])[0]
    dt = 1e-2
    times = np.array([i*dt for i in range(N)])
    reshaped_NW_33 = np.reshape(loaded_NW_33['lya'], (N, 256, 64))

    loaded_SW_33 = np.load('/home/x-loconnor/mhd/d2/quadrants/lya_SW_33_0033.npz')
    reshaped_SW_33 = np.reshape(loaded_SW_33['lya'], (N, 256, 64))

    loaded_SE_33 = np.load('/home/x-loconnor/mhd/d2/quadrants/lya_SE_33_0033.npz')
    reshaped_SE_33 = np.reshape(loaded_SE_33['lya'], (N, 256, 64))

    loaded_NE_33 = np.load('/home/x-loconnor/mhd/d2/quadrants/lya_NE_33_0033.npz')
    reshaped_NE_33 = np.reshape(loaded_NE_33['lya'], (N, 256, 64))

    reshaped = reshaped_SW_33
    offset_row, offset_col = 0, 0
    for row in range(256):
        for col in range(64):
            data_vec = reshaped[:, row, col]
            x, y = times, data_vec
            slope, intercept, r, p, std_err = stats.linregress(x, y)
            ftle2[row + offset_row, col + offset_col] = slope

    reshaped = reshaped_SE_33
    offset_row, offset_col = 0, 64
    for row in range(256):
        for col in range(64):
            data_vec = reshaped[:, row, col]
            x, y = times, data_vec
            slope, intercept, r, p, std_err = stats.linregress(x, y)
            ftle2[row + offset_row, col + offset_col] = slope

    reshaped = reshaped_NW_33
    offset_row, offset_col = 256, 0
    for row in range(256):
        for col in range(64):
            data_vec = reshaped[:, row, col]
            x, y = times, data_vec
            slope, intercept, r, p, std_err = stats.linregress(x, y)
            ftle2[row + offset_row, col + offset_col] = slope

    reshaped = reshaped_NE_33
    offset_row, offset_col = 256, 64
    for row in range(256):
        for col in range(64):
            data_vec = reshaped[:, row, col]
            x, y = times, data_vec
            slope, intercept, r, p, std_err = stats.linregress(x, y)
            ftle2[row + offset_row, col + offset_col] = slope


vmin = min(ftle1.min(), ftle2.min())
vmax = max(ftle1.max(), ftle2.max())

Np = 128
zn = np.linspace( 0, 2*np.pi, 4*Np+1, endpoint=False)[:-1] + 2*np.pi/(4*Np+1)
xn = np.linspace(-1,       1,   Np+1, endpoint=False)[:-1] + 2/(Np + 1)
positionsx = np.array([[xn[j] for j in range(Np)] for i in range(4*Np)])
positionsz = np.array([[zn[i] for j in range(Np)] for i in range(4*Np)])

plt.rcParams.update({'font.size': 16})
fig = plt.figure(figsize=(6, 7.5))
# fig = plt.figure()
# gs = GridSpec(1, 2, width_ratios=[1, 1], wspace=0, hspace=0,
#             left=0.1, right=0.9,
#             bottom=0.1, top=0.9)
# gs = GridSpec(1, 2, figure=fig, wspace=0, hspace=0)
# plot_positions = [0, 1]

gs = GridSpec(1, 3, width_ratios=[1, 1, 0.05], wspace=0)

# # Set same limits for both plots
# mesh1.set_clim(vmin, vmax)
# mesh2.set_clim(vmin, vmax)


ax = fig.add_subplot(gs[0, 0], adjustable='box', aspect=1)
pc = ax.pcolormesh(positionsx, positionsz, ftle1, cmap='inferno', rasterized=True, vmin=vmin, vmax=vmax)
ax.set_yticks([0, np.pi, 2*np.pi])
ax.set_yticklabels(["0", r"$\pi$", r"$2\pi$"])
ax.set_ylabel('z')
ax.set_xticks([-1, 0, 1])
ax.set_xlabel('x')
ax.set_title(r"$t=0$")
ax = fig.add_subplot(gs[0, 1], adjustable='box', aspect=1)
pc = ax.pcolormesh(positionsx, positionsz, ftle2, cmap='inferno', rasterized=True, vmin=vmin, vmax=vmax)
ax.set_xticks([-1, 0, 1])
ax.set_xlabel('x')
ax.set_yticks([0, np.pi, 2*np.pi])
# ax.set_ylabel('z')
ax.set_yticklabels([" ", " ", " "])
ax.set_title(r"$t=T/4$")
# cbar = fig.colorbar(im1, ax=[ax1, ax2], location='right', pad=0.02)
cax = fig.add_subplot(gs[2])
cbar = fig.colorbar(pc, cax=cax)
from matplotlib.ticker import FormatStrFormatter
cbar.formatter = FormatStrFormatter('%.2f')  # 2 decimal places
cbar.update_normal(pc)  # Update with the new formatter

plt.tight_layout()
plt.suptitle('FTLE')
fig.subplots_adjust(top=0.87)

dpi = 2400

figname = path + '/ftle.pdf'
plt.savefig(figname, dpi=dpi)
print(figname)

figname = path + '/ftle.png'
plt.savefig(figname, dpi=dpi)
print(figname)