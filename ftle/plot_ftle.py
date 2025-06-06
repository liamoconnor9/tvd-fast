import matplotlib.pyplot as plt
from scipy import stats
from dedalus import public as de
import numpy as np
import matplotlib.pyplot as plt
import sys

# loaded = np.load('/home/x-loconnor/ABC/FTLE_non_linear/lya.npz')
n = 1
loaded = np.load('/home/x-loconnor/mhd/d2/data30/lya_{}.npz'.format(str(n).zfill(4)))
# print(loaded['lya'].shape)
# print(np.max(loaded['lya'], axis=0).shape)
N = np.shape(loaded['lya'])[0]
dt = 1e-2
times = np.array([i*dt for i in range(N)])
reshaped = np.reshape(loaded['lya'], (N, 256, 64))
print(np.shape(reshaped))
data = np.zeros((256, 64))
counter = 1
for row in range(256):
    for col in range(64):
        counter += 1
        data_vec = reshaped[:, row, col]
        x, y = times, data_vec

        slope, intercept, r, p, std_err = stats.linregress(x, y)
        if counter % 1000 == 0:
            plt.plot(x, y)
            plt.show()
        data[row, col] = slope
        # print(slope)
        # print(std_err)

        # sys.exit()

# data_lin = np.max(loaded['lya'], axis=0)
# data = np.reshape(data_lin, (256, 64))

Lx = 2
Lz = 2*np.pi
x_basis = de.Chebyshev('x', 64, interval=(-Lx/2, Lx/2), dealias=3/2)
z_basis = de.Fourier('z',256, interval=(0, Lz), dealias=3/2)
domain = de.Domain([z_basis, x_basis], grid_dtype=np.float64)

z = domain.grid(0)
x = domain.grid(1)

fig, ax = plt.subplots()
pc = plt.pcolormesh(x, z, data)
cbar = fig.colorbar(pc)
ax.set_aspect(1)
plt.tight_layout()
plt.savefig('/home/x-loconnor/mhd/d2/data/figs/ftle30_{}.png'.format(str(n).zfill(4)))