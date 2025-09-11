import numpy as np
import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
import h5py
import dedalus.public as d3
from mpi4py import MPI
CW = MPI.COMM_WORLD
sys.path.append("..") # Adds higher directory to python modules path.
from vp_bvp_func import *
from glob import glob
import logging
logger = logging.getLogger(__name__)
from docopt import docopt
from pathlib import Path
from configparser import ConfigParser


if len(sys.argv) > 1:
    suffix = sys.argv[1]
else:
    logger.info('provide suffix')
    raise
filename = "{}/{}/options.cfg".format(path, suffix)

config = ConfigParser()
config.read(str(filename))
logger.info(config.items('parameters'))

suffix = eval(config.get('parameters', 'suffix'))
logger_cadence = config.getint('parameters','logger_cadence')
load_cp = eval(config.get('parameters','load_cp')) #what if we accepted string initial conditions if they start and end with $
isHydro = config.getboolean('parameters','isHydro')
isKinematic = config.getboolean('parameters','isKinematic')
doLoadTimestep = config.getboolean('parameters','doLoadTimestep')
Ny = config.getint('parameters','Ny')
Ly = eval(config.get('parameters','Ly'))
Nz = config.getint('parameters','Nz')
Lz = eval(config.get('parameters','Lz'))
Nx = config.getint('parameters','Nx')
Lx = eval(config.get('parameters','Lx'))
ic_scale_u = config.getfloat('parameters','ic_scale_u')
ic_scale_A = config.getfloat('parameters','ic_scale_A')
Ro = config.getfloat('parameters','Ro')
Re = config.getfloat('parameters','Re')
Rm = config.getfloat('parameters','Rm')
B0_coeff = config.getfloat('parameters', 'B0_coeff')
init_timestep = config.getfloat('parameters', 'init_timestep')
max_timestep = config.getfloat('parameters', 'max_timestep')
cfl_safety = config.getfloat('parameters', 'cfl_safety')
stop_sim_time = config.getfloat('parameters', 'stop_sim_time') + max_timestep
wall_time = 60. * 60. * config.getfloat('parameters', 'wall_time_hr')
timestepper = eval(config.get('parameters', 'timestepper'))
scalars_sim_dt = config.getfloat('parameters','scalars_sim_dt')
cp_sim_dt = config.getfloat('parameters','cp_sim_dt')
cp_scale = config.getfloat('parameters','cp_scale')
sp_sim_dt = config.getfloat('parameters','sp_sim_dt')
sp_scale = config.getfloat('parameters','sp_scale')

ary = Ly / Lx
arz = Lz / Lx
vol = Ly * Lz * Lx
Pm = Rm / Re
S = -1
B0_z = B0_coeff * (-Ro * Lx**2 * S * np.pi**(-2))

ncpu = MPI.COMM_WORLD.size
log2 = np.log2(ncpu)
if log2 == int(log2):
    mesh = [int(2**np.ceil(log2/2)),int(2**np.floor(log2/2))]
else:
    logger.error("pretty sure this shouldn't happen... log2(ncpu) is not an int?")
    
logger.info("running on processor mesh={}".format(mesh))

# Bases
coords = d3.CartesianCoordinates('y', 'z', 'x')
dealias = 3/2
dist = d3.Distributor(coords, dtype=np.float64, mesh=mesh)
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
integy = lambda A: d3.Integrate(A, 'y')
integz = lambda A: d3.Integrate(A, 'z')
integx = lambda A: d3.Integrate(A, 'x')
integ = lambda A: integy(integz(integx(A)))
dy = lambda A: d3.Differentiate(A, coords['y'])
dz = lambda A: d3.Differentiate(A, coords['z'])
dx = lambda A: d3.Differentiate(A, coords['x'])
from dedalus.core.operators import TimeDerivative
dt = lambda argy: TimeDerivative(argy)
lift_basis = xbasis.clone_with(a=1/2, b=1/2) # First derivative basis
lift = lambda A, n: d3.Lift(A, lift_basis, n)

# nccs
U0 = dist.VectorField(coords, name='U0', bases=xbasis)
U0['g'][0] = S * x
z_hat = dist.VectorField(coords, name='z_hat', bases=xbasis)
z_hat['g'][1] = 1

t = dist.Field(name='t')
p = dist.Field(name='p', bases=bases)
u = dist.VectorField(coords, name='u', bases=bases)
uinit = dist.VectorField(coords, name='uinit', bases=bases)
taup = dist.Field(name='taup')
tau1u = dist.VectorField(coords, name='tau1u', bases=(ybasis,zbasis))
tau2u = dist.VectorField(coords, name='tau2u', bases=(ybasis,zbasis))
grad_u = d3.grad(u) + ex*lift(tau1u,-1)

A = dist.VectorField(coords, name='A', bases=bases)
Ainit = dist.VectorField(coords, name='Ainit', bases=bases)
phi = dist.Field(name='phi', bases=bases)
tau1A = dist.VectorField(coords, name='tau1A', bases=(ybasis,zbasis))
tau2A = dist.VectorField(coords, name='tau2A', bases=(ybasis,zbasis))
grad_A = d3.grad(A) + ex*lift(tau1A,-1)
b = d3.Curl(A)

# Initial conditions
pert_scales = 0.25
fh_mode = 'overwrite'
imported_time = 0.0

load_cp = suffix + '/checkpoint/checkpoint_s21.h5'
load_path = "{}/{}".format(path, load_cp)
with h5py.File(load_path, "r") as file:
    u_temp = ic_scale_u * file['tasks']['u'][()][0, ...]
    u.load_from_global_grid_data(u_temp)
    u.change_scales(1)
    A.load_from_global_grid_data(ic_scale_A * file['tasks']['A'][()][0, ...])
    A.change_scales(1)
                # Ainit['g'] = A['g'].copy()

A0 = A.copy()

load_cp = suffix + '/checkpoint/checkpoint_s22.h5'
load_path = "{}/{}".format(path, load_cp)
with h5py.File(load_path, "r") as file:
    u_temp = ic_scale_u * file['tasks']['u'][()][0, ...]
    u.load_from_global_grid_data(u_temp)
    u.change_scales(1)
    A.load_from_global_grid_data(ic_scale_A * file['tasks']['A'][()][0, ...])
    A.change_scales(1)
                # Ainit['g'] = A['g'].copy()

A1 = A.copy()

def inner_p(field1, field2):
    return integ(field1 @ field2)

def get_norm(field):
    return np.sqrt(inner_p(field, field))


norm_A0 = get_norm(A0).evaluate()
A0_normed = A0 / norm_A0
norm_A1 = get_norm(A1).evaluate()
A1_normed = A1 / norm_A1

proj = inner_p(A0_normed, A1_normed).evaluate()


# ratio = (A1 / A0).evaluate()
if CW.rank == 0:
    ratio = norm_A1['g'] / norm_A0['g']
    lamb = np.log(ratio) / 23.33
    print("normalized projection (ideally = 1): {}".format(proj['g']))
    print("norm_A0: {}".format(norm_A0['g']))
    print("norm_A1: {}".format(norm_A1['g']))
    print("growth multiplier (ratio of norms): {}".format(ratio))
    print("lamb (growth rate): {}".format(lamb))
