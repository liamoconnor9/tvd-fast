"""
Run rotating plane couette flow hydro sim and calculate FTLEs

Usage:
    shear_paths.py [options]

Options:
    --restart_file=<restart_file>     Location of restart file [default: None]
    --restartN=<restartN>             Restart number [default: -1]
    --T=<T>                           Run time [default: 200.0]
    --N=<N>                           Number of particles in each direction [default: 64]
    --mesh=<mesh>                     Processor mesh for 3-D runs
    --name=<name>                     run name. specify quadrant (e.g., NE in name) [default: None]
"""

import numpy as np
import time
import h5py
from mpi4py import MPI
import copy
from dedalus import public as de
from dedalus.extras import flow_tools
import os
from tools import load_state_partial
import logging
logger = logging.getLogger(__name__)
from docopt import docopt
import sys
sys.path.append('../particle_tracker/')
import particles as particles
import os
path = os.path.dirname(os.path.abspath(__file__))

# parse arguments
args = docopt(__doc__)

# Parameters
Reynolds     = 150
Ro           = 3
restartN     = int(args['--restartN'])
T       = float(args['--T'])
N       = int(args['--N'])
name       = args['--name']
if name == 'None':
    print('specify name')
    raise

restart_file = path + "/../kin_2p6d_Rm1500_RSTRT5/checkpoint/checkpoint_s{}.h5".format(restartN)

logger.info("####################################################################")
logger.info("restartN = {}".format(restartN))
logger.info("####################################################################")

comm = MPI.COMM_WORLD
rank = comm.rank
ncpu = comm.size
mesh = args['--mesh']
if mesh is not None:
    mesh = mesh.split(',')
    mesh = [int(mesh[0]), int(mesh[1])]
    logger.info("running on processor mesh={}".format(mesh))
else:
    log2 = np.log2(ncpu)
    if log2 == int(log2):
        mesh = [int(2**np.ceil(log2/2)),int(2**np.floor(log2/2))]
    logger.info("running on processor mesh={}".format(mesh))

logger.info('Parameters: Re: %e, Ro: %e' %(Reynolds, Ro))
logger.info('Run time  : T:  %e' %(T))

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#Aspect ratio
Lx, Ly, Lz = (2, 2*np.pi, 2*np.pi)
nx, ny, nz = (64, 16, 256)

# Create bases and domain
x_basis = de.Chebyshev('x', nx, interval=(-Lx/2, Lx/2), dealias=3/2)
y_basis = de.Fourier('y',ny, interval=(0, Ly), dealias=3/2)
z_basis = de.Fourier('z',nz, interval=(0, Lz), dealias=3/2)
domain = de.Domain([y_basis, z_basis, x_basis], grid_dtype=np.float64,mesh=mesh)
slices = domain.distributor.grid_layout.slices(scales=1)

logger.info('Purely hydro simulation')
problem = de.IVP(domain, variables=['u','v','w','ux','vx','wx','p'])

problem.parameters['ReInv'] = 1./Reynolds
problem.parameters['RoInv'] = 1./Ro
problem.parameters['Ro'] = Ro

problem.substitutions['Lap(A)'] = "dx(dx(A)) + dy(dy(A)) + dz(dz(A))"
problem.substitutions['ad(A)'] = "u*dx(A) + v*dy(A) + w*dz(A)"

# cross(zhat, {u, v, w}) = {-v, u, 0}

# Hydrodynamic
problem.add_equation("ux + dy(v) + dz(w) = 0")
problem.add_equation("dt(u) - ReInv*(dx(ux) + dy(dy(u)) + dz(dz(u))) +  dx(p) + RoInv*(-v)   = - u*ux - v*dy(u) - w*dz(u)")
problem.add_equation("dt(v) - ReInv*(dx(vx) + dy(dy(v)) + dz(dz(v))) +  dy(p) + RoInv*( u)   = - u*vx - v*dy(v) - w*dz(v)")
problem.add_equation("dt(w) - ReInv*(dx(wx) + dy(dy(w)) + dz(dz(w))) +  dz(p)                = - u*wx - v*dy(w) - w*dz(w)")
problem.add_equation("ux - dx(u) = 0")
problem.add_equation("vx - dx(v) = 0")
problem.add_equation("wx - dx(w) = 0")

problem.add_bc("left(u) = 0")
problem.add_bc("left(v) = -left(x)")
problem.add_bc("left(w) = 0")
problem.add_bc("right(u) = 0", condition="(ny != 0) or (nz != 0)")
problem.add_bc("right(v) = -right(x)")
problem.add_bc("right(w) = 0")
problem.add_bc("integ_x(p) = 0", condition="(ny == 0) and (nz == 0)")

ts = de.timesteppers.RK443
solver =  problem.build_solver(ts)

y = domain.grid(0)
z = domain.grid(1)
x = domain.grid(2)

u = solver.state['u']
v = solver.state['v']
w = solver.state['w']
ux = solver.state['ux']
vx = solver.state['vx']
wx = solver.state['wx']

if(restart_file != 'None'):
    logger.info('Restarting from: %s' % restart_file)
    try:
        with h5py.File(restart_file, "r") as f:
            uvec = f['tasks']['u'][0, :, 0, ...]
            uvec_y = uvec[0, ...]
            uvec_z = uvec[1, ...]
            uvec_x = uvec[2, ...]
            v['g'] = uvec_y[slices[1:]]
            w['g'] = uvec_z[slices[1:]]
            u['g'] = uvec_x[slices[1:]]


    except:
        logger.info('Run non-linear simmulation to get restart states first')
        raise
    T += solver.sim_time

try:
    write, last_dt = load_state_partial(solver,restart_file, 0)
except:
    logger.info('Wrong file name. Make sure you run a non-linear simmulation to get restart states first')
    raise

solver.sim_time = 0.0
solver.stop_sim_time = T
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

initial_dt = 1e-3
cfl = flow_tools.CFL(solver,initial_dt,safety=0.8)
cfl.add_velocities(('u','v','w'))

flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("sqrt(u*u + v*v + w*w) / ReInv", name='Re')

logger.info('Starting loop')
start_time = time.time()

# Initiate particles
particles = particles.particles(4*N**2,domain)

# Rewrite initial positions with equispaced locations
# n = int(np.sqrt(particles.N))
if "NW" in name:
    zmin = np.pi
    zmax = 2*np.pi
    xmin = -1
    xmax = 0
elif "NE" in name:
    zmin = np.pi
    zmax = 2*np.pi
    xmin = 0
    xmax = 1
elif "SW" in name:
    zmin = 0
    zmax = np.pi
    xmin = -1
    xmax = 0
elif "SE" in name:
    zmin = 0
    zmax = np.pi
    xmin = 0
    xmax = 1

zn = np.linspace( zmin, zmax, 4*N+1)[:-1]
xn = np.linspace( xmin, xmax,   N+1)[:-1]

particles.positions = np.array([(np.pi, zn[i],xn[j]) for i in range(4*N) for j in range(N)])

# To Save memory must perform some of the FTLE calculation here
lya = []
times = []

while solver.ok:
    dt = 0.01
    solver.step(dt)

    particles.step(dt,(v,w,u))
    particles.stepStress(dt,(v,w,u))
    if solver.iteration % 1 == 0:
        # Lyapunov exponent calculation
        C = np.einsum('kji,kij->kij',particles.J,particles.J)
        lya.append(np.log(np.einsum('kii->k',C)/2.)/2.)
        times.append(solver.sim_time)

    if (solver.iteration-1) % 10 == 0:
        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))

end_time = time.time()

# Save data
folder_name = 'quadrants'

if MPI.COMM_WORLD.rank == 0:
    if not os.path.exists('{}/'.format(folder_name)):
        os.mkdir(folder_name)

if(rank==0):
    np.savez('{}/lya_{}_{}'.format(folder_name, name, str(restartN).zfill(4)), lya=np.array(lya),times=times, positions=particles.positions)
logger.info('shape lya = {}'.format(np.shape(np.array(lya))))

# Print statistics
logger.info('Run time: %f' %(end_time-start_time))
logger.info('Iterations: %i' %solver.iteration)
