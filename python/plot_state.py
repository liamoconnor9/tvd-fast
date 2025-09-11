from scipy import stats
from dedalus import public as d3
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
import os
path = os.path.dirname(os.path.abspath(__file__))
from matplotlib.gridspec import GridSpec
import h5py

half = False
if len(sys.argv) == 2:
    N1 = int(sys.argv[1])
    auto = True
else:
    N1 = 1
    auto = False
N2 = (N1 + 32) % 128

if N1 == 0 or N2 == 0:
    sys.exit()

suffix = "kin_2p6d_Rm1500_RSTRT5"
Nz, Nx = 256, 64

with h5py.File("{}/{}/checkpoint/checkpoint_s{}.h5".format(path, suffix, N1), "r") as f:
    uvec = f['tasks']['u'][0, :, 0, ...]
    sim_time = f['scales']['sim_time'][0]
    uvec_y = uvec[0, ...]
    uvec_z = uvec[1, ...]
    uvec_x = uvec[2, ...]

    bvec = f['tasks']['b'][0, :, 0, ...]
    bvec_y = bvec[0, ...]
    bvec_z = bvec[1, ...]
    bvec_x = bvec[2, ...]

    vy_t1 = uvec_y
    vz_t1 = uvec_z
    vx_t1 = uvec_x

    by_t1_e3 = bvec_y
    bz_t1 = bvec_z
    bx_t1 = bvec_x
    

with h5py.File("{}/{}/checkpoint/checkpoint_s{}.h5".format(path, suffix, N2), "r") as f:
    uvec = f['tasks']['u'][0, :, 0, ...]
    uvec_y = uvec[0, ...]
    uvec_z = uvec[1, ...]
    uvec_x = uvec[2, ...]

    bvec = f['tasks']['b'][0, :, 0, ...]
    bvec_y = bvec[0, ...]
    bvec_z = bvec[1, ...]
    bvec_x = bvec[2, ...]

    vy_t2 = uvec_y
    vz_t2 = uvec_z
    vx_t2 = uvec_x

    by_t2_e3 = bvec_y
    bz_t2 = bvec_z
    bx_t2 = bvec_x


Lx = 2
Lz = 2*np.pi
coords = d3.CartesianCoordinates('z', 'x')
dealias = 3/2
dist = d3.Distributor(coords, dtype=np.float64)
zbasis = d3.RealFourier(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)
xbasis = d3.ChebyshevT(coords['x'], size=Nx, bounds=(-Lx / 2.0, Lx / 2.0), dealias=dealias)
z = dist.local_grid(zbasis)
x = dist.local_grid(xbasis)
X_e3, Y_e3 = np.meshgrid(x, z)

temp = dist.Field(name='temp', bases=[zbasis, xbasis])
temp['g'] = by_t1_e3
temp.change_scales(8)
by_t1_e3 = temp['g'].copy()
temp.change_scales(1)

temp['g'] = by_t2_e3
temp.change_scales(8)
by_t2_e3 = temp['g'].copy()
temp.change_scales(1)

suffix = "Ly22_kin_Rm150000_RSTRT1"
Nz, Nx = 2048, 512

with h5py.File("{}/{}/checkpoint/checkpoint_s{}.h5".format(path, suffix, N1), "r") as f:
    uvec = f['tasks']['u'][0, :, 0, ...]
    sim_time = f['scales']['sim_time'][0]
    uvec_y = uvec[0, ...]
    uvec_z = uvec[1, ...]
    uvec_x = uvec[2, ...]

    bvec = f['tasks']['b'][0, :, 0, ...]
    bvec_y = bvec[0, ...]
    bvec_z = bvec[1, ...]
    bvec_x = bvec[2, ...]

    vy_t1 = uvec_y
    vz_t1 = uvec_z
    vx_t1 = uvec_x

    by_t1_e5 = bvec_y
    bz_t1 = bvec_z
    bx_t1 = bvec_x
    

with h5py.File("{}/{}/checkpoint/checkpoint_s{}.h5".format(path, suffix, N2), "r") as f:
    uvec = f['tasks']['u'][0, :, 0, ...]
    uvec_y = uvec[0, ...]
    uvec_z = uvec[1, ...]
    uvec_x = uvec[2, ...]

    bvec = f['tasks']['b'][0, :, 0, ...]
    bvec_y = bvec[0, ...]
    bvec_z = bvec[1, ...]
    bvec_x = bvec[2, ...]

    vy_t2 = uvec_y
    vz_t2 = uvec_z
    vx_t2 = uvec_x

    by_t2_e5 = bvec_y
    bz_t2 = bvec_z
    bx_t2 = bvec_x


Lx = 2
Lz = 2*np.pi
coords = d3.CartesianCoordinates('z', 'x')
dealias = 3/2
dist = d3.Distributor(coords, dtype=np.float64)
zbasis = d3.RealFourier(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)
xbasis = d3.ChebyshevT(coords['x'], size=Nx, bounds=(-Lx / 2.0, Lx / 2.0), dealias=dealias)
z = dist.local_grid(zbasis)
x = dist.local_grid(xbasis)
X_e5, Y_e5 = np.meshgrid(x, z)


ftle_source = 'data_thirtyone'
# ftle_source = 'data30'
loaded = np.load('/home/x-loconnor/mhd/d2/{}/lya_{}.npz'.format(ftle_source, str(N1).zfill(4)))
N = np.shape(loaded['lya'])[0]
dt = 1e-2
times = np.array([i*dt for i in range(N)])
reshaped = np.reshape(loaded['lya'], (N, 256, 64))

Np = 64
zn = np.linspace( 0, 2*np.pi, 4*Np+1)[:-1]
xn = np.linspace(-1,       1,   Np+1)[:-1]
positionsx = np.array([[xn[j] for j in range(Np)] for i in range(4*Np)])
positionsz = np.array([[zn[i] for j in range(Np)] for i in range(4*Np)])


ftle1 = np.zeros((256, 64))
for row in range(256):
    for col in range(64):
        data_vec = reshaped[:, row, col]
        x, y = times, data_vec
        slope, intercept, r, p, std_err = stats.linregress(x, y)
        ftle1[row, col] = slope

loaded = np.load('/home/x-loconnor/mhd/d2/{}/lya_{}.npz'.format(ftle_source, str(N2).zfill(4)))
reshaped = np.reshape(loaded['lya'], (N, 256, 64))
ftle2 = np.zeros((256, 64))
for row in range(256):
    for col in range(64):
        data_vec = reshaped[:, row, col]
        x, y = times, data_vec
        slope, intercept, r, p, std_err = stats.linregress(x, y)
        ftle2[row, col] = slope

by_t1_e3 = by_t1_e3 / np.max(by_t1_e3)
by_t2_e3 = by_t2_e3 / np.max(by_t2_e3)

by_t1_e5 = by_t1_e5 / np.max(by_t1_e5)
by_t2_e5 = by_t2_e5 / np.max(by_t2_e5)

# Six different patterns
patterns = [
    vy_t1,
    by_t1_e3,
    by_t1_e5,
    vy_t2,
    by_t2_e3,
    by_t2_e5,
]

cmaps = ['RdBu_r', 'PiYG', 'PiYG', 'RdBu_r', 'PiYG', 'PiYG']

# Create figure with tight layout

plt.rcParams.update({'font.size': 9})
# plt.figure(figsize=(6, 4))


# GridSpec with no spacing between subplots
if not half:
    fig = plt.figure()
    gs = GridSpec(1, 7, width_ratios=[1, 1, 1, 0.2, 1, 1, 1], wspace=0, hspace=0,
                left=0.1, right=0.95,
                bottom=0.1, top=0.9)
    plot_positions = [0, 1, 2, 4, 5, 6]
else:
    fig = plt.figure(figsize=(7, 7))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0, hspace=0,
                left=0.1, right=0.9,
                bottom=0.1, top=0.9)
    plot_positions = [0, 1, 2]

# Create subplots
pcs = []
# titles = [r"$\mathbf{u} \cdot \mathbf{\hat{y}}$", r"$\mathbf{b} \cdot \mathbf{\hat{y}}$" + "\n" + r"$\rm{Rm}=1.5\cdot 10^3$", r"$\mathbf{b} \cdot \mathbf{\hat{y}}$" + "\n" + r"$\rm{Rm}=1.5\cdot 10^5$", r"$\mathbf{u} \cdot \mathbf{\hat{y}}$", r"$\mathbf{b} \cdot \mathbf{\hat{y}}$" + "\n" + r"$\rm{Rm}=1.5\cdot 10^3$", r"$\mathbf{b} \cdot \mathbf{\hat{y}}$" + "\n" + r"$\rm{Rm}=1.5\cdot 10^5$"]
titles = [r"$\mathbf{u} \cdot \mathbf{\hat{y}}$", r"$\mathbf{b} \cdot \mathbf{\hat{y}}$", r"$\mathbf{b} \cdot \mathbf{\hat{y}}$", r"$\mathbf{u} \cdot \mathbf{\hat{y}}$", r"$\mathbf{b} \cdot \mathbf{\hat{y}}$", r"$\mathbf{b} \cdot \mathbf{\hat{y}}$"]


# Create subplots
annotations = ["A", "B", "C", "D", "E", "F"]
for i, pos in enumerate(plot_positions):
    ax = fig.add_subplot(gs[pos], adjustable='box', aspect=1)
    # if i == 1 or i == 4:
    #     xplt, yplt = X_e3, Y_e3
    # else:
    #     xplt, yplt = X_e5, Y_e5
    xplt, yplt = X_e5, Y_e5

    pc = ax.pcolormesh(xplt, yplt, patterns[i], cmap=cmaps[i], rasterized=True)
    # ax.text(annotations[i], [1, 1])
    ax.text(0.45, 0.93, annotations[i], bbox=dict(facecolor='white', alpha=1.0), transform=ax.transAxes)
    pcs.append(pc)

    ax.set_xticks([])
    ax.set_yticks([0, np.pi, 2*np.pi])
    if i == 0:
        ax.set_yticklabels(["0", r"$\pi$", r"$2\pi$"])
        ax.set_ylabel('z')
    else:
        ax.set_yticklabels([" ", " ", " "])

    if i == 0 or i == 3:
        ax.set_xticks([-1, 0])
        ax.set_xticklabels(["-1", "0"])    
        ax.set_xlabel('x')
    elif i == 1 or i == 4:
        ax.set_xticks([-1, 0])
        ax.set_xticklabels(["-1", "0"])    
        ax.set_xlabel('x')
    else:
        ax.set_xticks([-1, 0, 1])
        ax.set_xticklabels(['-1', "0", '1'])    
        ax.set_xlabel('x')
    if half:
        ax.set_xticks([-1, 0, 1])
        ax.set_xticklabels(['-1', "0", '1'])    
        ax.set_xlabel('x')

    ax.set_title(titles[i])

if not half:
    cby = 0.09
else:
    cby = 0.075

if not half:
    # Add colorbar
    cbar_ax = fig.add_axes([0.1, cby, 0.4, 0.015])  # [left, bottom, width, height]
    cbar = fig.colorbar(pcs[0], cax=cbar_ax, orientation='horizontal')
    # cbar.set_label('Value Scale')
    cbar_ax.xaxis.set_ticks_position('bottom')
    cbar_ax.set_xticks([-1, 0, 1])
    cbar_ax.set_xticklabels(["-1", "0", "1"])

    cbar_ax = fig.add_axes([0.55, cby, 0.4, 0.015])  # [left, bottom, width, height]
    cbar = fig.colorbar(pcs[2], cax=cbar_ax, orientation='horizontal')
    # cbar.set_label('Value Scale')
    cbar_ax.xaxis.set_ticks_position('bottom')
    cbar_ax.set_xticks([-1, 0, 1])
    cbar_ax.set_xticklabels(["-1", "0", "1"])

else:
    # Add colorbar
    cbx = 0.92
    cbar_ax = fig.add_axes([cbx, 0.7, 0.03, 0.1])  # [left, bottom, width, height]
    cbar = fig.colorbar(pcs[0], cax=cbar_ax, orientation='vertical')
    cbar_ax.set_yticks([-1, 0, 1])
    cbar_ax.set_yticklabels(["-1", "0", "1"])

    cbar_ax = fig.add_axes([cbx, 0.1, 0.03, 0.1])  # [left, bottom, width, height]
    cbar = fig.colorbar(pcs[2], cax=cbar_ax, orientation='vertical')
    cbar_ax.set_yticks([-1, 0, 1])
    cbar_ax.set_yticklabels(["-1", "0", "1"])
    plt.suptitle('t = {:.4f}'.format(sim_time))

if not half:
    xcoord = -3.7
    ycoord = 1.15
    yshift = 0.1
    ax.annotate(r"$t=0$", xy=(xcoord, ycoord), xytext=(xcoord, ycoord + yshift), xycoords='axes fraction', 
                ha='center', va='top',
                bbox=dict(boxstyle='square', fc='white', color='k'),
                arrowprops=dict(arrowstyle='-[, widthB=9.0, lengthB=0.5', lw=2.0, color='k'))

    xcoord = -0.5
    # ycoord = -0.2
    ax.annotate(r"$t=T/4$", xy=(xcoord, ycoord), xytext=(xcoord, ycoord + yshift), xycoords='axes fraction', 
                ha='center', va='top',
                bbox=dict(boxstyle='square', fc='white', color='k'),
                arrowprops=dict(arrowstyle='-[, widthB=9.0, lengthB=0.5', lw=2.0, color='k'))

# plt.tight_layout()
filetype = '.pdf'
dpi = 2400
if half:
    if auto:
        plt.savefig('{}/halves/half_{}{}'.format(path, str(N1).zfill(4)), filetype, dpi=dpi)
    else:
        plt.savefig('{}/half{}'.format(path, filetype), dpi=dpi)
else:
    if auto:
        plt.savefig('{}/states/state_{}{}'.format(path, str(N1).zfill(4), filetype), dpi=dpi)
    else:
        plt.savefig('{}/state{}'.format(path, filetype), dpi=dpi)
print('done state.pdf')
filetype = '.png'
if half:
    if auto:
        plt.savefig('{}/halves/half_{}{}'.format(path, str(N1).zfill(4)), filetype, dpi=dpi)
    else:
        plt.savefig('{}/half{}'.format(path, filetype), dpi=dpi)
else:
    if auto:
        plt.savefig('{}/states/state_{}{}'.format(path, str(N1).zfill(4), filetype), dpi=dpi)
    else:
        plt.savefig('{}/state{}'.format(path, filetype), dpi=dpi)

print('done state.png')