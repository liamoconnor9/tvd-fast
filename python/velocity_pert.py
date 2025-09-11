import numpy as np
import dedalus.public as d3
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from mpi4py import MPI
CW = MPI.COMM_WORLD
import os
path = os.path.dirname(os.path.abspath(__file__))
from docopt import docopt
from configparser import ConfigParser
from pathlib import Path
import logging
logger = logging.getLogger(__name__)


suffix = 'pm25_nu1en3'

try:    
    args = docopt(__doc__)
    filename = Path(args['<config_file>'])
    seed = int(args['<seed>'])
except:
    logger.warning('invalid config file supplied. using default... w/ seed = 1')
    filename = "{}/{}/mri_options.cfg".format(path, suffix)
    seed = 1

config = ConfigParser()
config.read(str(filename))

logger.info(config.items('parameters'))
scale = config.getfloat('parameters','scale')
ra_scale = config.getfloat('parameters','ra_scale')
suffix = eval(config.get('parameters', 'suffix'))
logger.info(suffix)
load_cp = config.getboolean('parameters','load_cp')
solve_EVP = config.getboolean('parameters','solve_EVP')
logger_cadence = config.getint('parameters','logger_cadence')
Ny = config.getint('parameters','Ny')
Ly = eval(config.get('parameters','Ly'))
Nz = config.getint('parameters','Nz')
Lz = eval(config.get('parameters','Lz'))
Nx = config.getint('parameters','Nx')
Lx = eval(config.get('parameters','Lx'))
tau = config.getfloat('parameters','tau')
nu = config.getfloat('parameters','nu')
Pm = config.getfloat('parameters','Pm')
eta = nu / Pm
Ro = config.getfloat('parameters', 'Ro')
f = config.getfloat('parameters', 'f')
S = -f * Ro
B0_coeff = config.getfloat('parameters', 'B0_coeff')
B0_z = B0_coeff * (-f * Lx**2 * S * np.pi**(-2))
logger.info('[INPUTTED] Ro = {}'.format(Ro))
logger.info('[INPUTTED] B0_coeff = {}'.format(B0_coeff))
logger.info('[OUTPUTTED] S = {}'.format(S))
logger.info('[OUTPUTTED] B0_z = {}'.format(B0_z))

# Bases
dtype = np.float64
coords = d3.CartesianCoordinates('y', 'z', 'x')
dist = d3.Distributor(coords, dtype=dtype)
dealias = 3/2

ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(0, Ly), dealias=dealias)
zbasis = d3.RealFourier(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)
xbasis = d3.ChebyshevT(coords['x'], size=Nx, bounds=(-Lx / 2.0, Lx / 2.0), dealias=dealias)

bases = (ybasis,zbasis,xbasis)
y = dist.local_grid(ybasis)
z = dist.local_grid(zbasis)
x = dist.local_grid(xbasis)
ey = dist.VectorField(coords, name='ey')
ez = dist.VectorField(coords, name='ez')
ex = dist.VectorField(coords, name='ex')
ey['g'][0] = 1
ez['g'][1] = 1
ex['g'][2] = 1


fz_hat = dist.VectorField(coords, name='fz_hat', bases=xbasis)
fz_hat['g'][1] = f

# damping timescale
TAU = dist.Field(name='TAU')

# Fields
p = dist.Field(name='p', bases=bases)
phi = dist.Field(name='phi', bases=bases)
u = dist.VectorField(coords, name='u', bases=bases)
psi = dist.VectorField(coords, name='psi', bases=bases)
binit = dist.VectorField(coords, name='binit', bases=bases)
A = dist.VectorField(coords, name='A', bases=bases)
taup = dist.Field(name='taup')
tau1u = dist.VectorField(coords, name='tau1u', bases=(ybasis,zbasis))
tau2u = dist.VectorField(coords, name='tau2u', bases=(ybasis,zbasis))
tau1A = dist.VectorField(coords, name='tau1A', bases=(ybasis,zbasis))
tau2A = dist.VectorField(coords, name='tau2A', bases=(ybasis,zbasis))

target_dir = "{}/{}/data/1/checkpoint/*.h5".format(path, suffix)
files = glob.glob(target_dir)

# nccs
U0 = dist.VectorField(coords, name='U0', bases=xbasis)
U0['g'][0] = S * x

def energy3D(arg):
    return d3.Integrate(d3.Integrate(d3.Integrate(arg @ arg / 2 / Ly / Lz / Lx, 'y'), 'z'), 'x').evaluate()['g'][0]

time_vec = []
ufull_evec = []
ubar_evec  = []
uprime_evec = []

for file in files:
    print(file)
    with h5py.File(file, "r") as f:
        time = f['scales']['sim_time'][()].squeeze()
        time_vec.append(time)
        udata = f['tasks']['u'][()][0, ...]
        u.change_scales(1)
        u['g'] = udata.copy()
        ufull = u + U0
        ubar = d3.Integrate(d3.Integrate(ufull, 'y'), 'z') / Ly / Lz
        ufull_energy = energy3D(ufull)
        ubar_energy = (d3.Integrate(ubar @ ubar / 2, 'x') / Lx).evaluate()['g'][0]
        uprime_energy = energy3D(ufull - ubar)

        ufull_evec.append(ufull_energy)
        ubar_evec.append(ubar_energy)
        uprime_evec.append(uprime_energy)


def get_title():
    return r"$\tau=$" + str(tau) + r"; $Pm=$" + str(Pm) + r"; $\nu=$" + str(nu)

plt.scatter(time_vec, ufull_evec, label = r"$||u||$", c='black')
plt.scatter(time_vec, ubar_evec, label = r"$||\langle u \rangle_{y, z}||$", marker='+', c='lime')
plt.scatter(time_vec, uprime_evec, label = r"$||u - \langle u \rangle_{y, z}||$", marker='x', c='purple')
plt.xlabel('time')
plt.title(get_title())
plt.legend()
figname = "{}/{}/decomp_lin.png".format(path, suffix)
print(figname)
plt.savefig(figname)
plt.close()

plt.scatter(time_vec, ufull_evec, label = r"$||u||$", c='black')
plt.scatter(time_vec, ubar_evec, label = r"$||\langle u \rangle_{y, z}||$", marker='+', c='lime')
plt.scatter(time_vec, uprime_evec, label = r"$||u - \langle u \rangle_{y, z}||$", marker='x', c='purple')
plt.xlabel('time')
plt.yscale('log')
plt.title(get_title())
plt.legend()
figname = "{}/{}/decomp_log.png".format(path, suffix)
print(figname)
plt.savefig(figname)