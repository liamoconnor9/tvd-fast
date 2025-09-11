"""
3D cartesian Rotating Plane Couette Flow (RPCF) Floquet analysis of the hydro system
Usage:
    floq_hydro.py <config_file>
"""

import numpy as np
import time

from dedalus import public as d3
from mpi4py import MPI
import logging
from docopt import docopt
logger = logging.getLogger(__name__)
from configparser import ConfigParser
from pathlib import Path
import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
import h5py
CW = MPI.COMM_WORLD
import logging
import pathlib
from glob import glob

try:    
    args = docopt(__doc__)
    filename = Path(args['<config_file>'])
except:
    filename = path + "/ecs5/mri_options.cfg"
config = ConfigParser()
config.read(str(filename))

logger.info(config.items('parameters'))

scale = config.getfloat('parameters','scale')
logger_cadence = config.getint('parameters','logger_cadence')

suffix = eval(config.get('parameters', 'suffix'))
load_cp = eval(config.get('parameters','load_cp'))
isHydro = config.getboolean('parameters','isHydro')
isKinematic = config.getboolean('parameters','isKinematic')
loadFromPert = config.getboolean('parameters','loadFromPert')
doLoadTimestep = config.getboolean('parameters','doLoadTimestep')

Ny = config.getint('parameters','Ny')
Ly = eval(config.get('parameters','Ly'))

Nz = config.getint('parameters','Nz')
Lz = eval(config.get('parameters','Lz'))

Nx = config.getint('parameters','Nx')
Lx = eval(config.get('parameters','Lx'))

ic_scale_u = config.getfloat('parameters','ic_scale_u')
ic_scale_A = config.getfloat('parameters','ic_scale_A')
nu = config.getfloat('parameters','nu')
Pm = config.getfloat('parameters','Pm')
Pr = config.getfloat('parameters','Pr')
tau = config.getfloat('parameters','tau')
eta = nu / Pm
evp_tol = config.getfloat('parameters','evp_tol')
NEV = config.getint('parameters','NEV')

isConvecting = config.getboolean('parameters', 'isConvecting')
isContinuation = config.getboolean('parameters', 'isContinuation')
if isContinuation:
    t_vec = eval(eval(config.get('parameters', 't_vec')))
    Ro_vec = eval(eval(config.get('parameters', 'Ro_vec')))
    Ra_vec = eval(eval(config.get('parameters', 'Ra_vec')))
    Ro = Ro_vec[0]
    Ra = Ra_vec[0]

isNoSlip = config.getboolean('parameters','isNoSlip')
f = config.getfloat('parameters', 'f')

S = -f * Ro
ary = Ly / Lx
arz = Lz / Lx
vol = Ly * Lz * Lx

# Evolution params
init_timestep = max_timestep = config.getfloat('parameters', 'timestep')
cfl_safety = config.getfloat('parameters', 'cfl_safety')
stop_sim_time = config.getfloat('parameters', 'stop_sim_time') + 0.01
wall_time = 60. * 60. * config.getfloat('parameters', 'wall_time_hr')
timestepper = eval(config.get('parameters', 'timestepper'))

# i/o
crumb_sim_dt = config.getfloat('parameters','crumb_sim_dt')
crumb_scale = config.getfloat('parameters','crumb_scale')

scalars_sim_dt = config.getfloat('parameters','scalars_sim_dt')
cp_sim_dt = config.getfloat('parameters','cp_sim_dt')
cp_scale = config.getfloat('parameters','cp_scale')
sp_sim_dt = config.getfloat('parameters','sp_sim_dt')

T = config.getfloat('parameters', 'T')
N = config.getint('parameters', 'N')

ncpu = CW.size
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

sim_index = 0
timestep = T / N

p = dist.Field(name='p', bases=bases)
u = dist.VectorField(coords, name='u', bases=bases)
u0 = dist.VectorField(coords, name='u0', bases=bases)
fz_hat = dist.VectorField(coords, name='fz_hat', bases=xbasis)
fz_hat['g'][1] = f

t = dist.Field()
tau_p = dist.Field()
tau1u = dist.VectorField(coords, name='tau1u', bases=(ybasis,zbasis))
tau2u = dist.VectorField(coords, name='tau2u', bases=(ybasis,zbasis))

lift_basis = xbasis.clone_with(a=1/2, b=1/2) # First derivative basis
lift = lambda A, n: d3.Lift(A, lift_basis, n)
grad_u = d3.grad(u) + ex*lift(tau1u,-1)

problem = d3.IVP([p, u, tau_p, tau1u, tau2u], time=t, namespace=locals())
problem.add_equation("trace(grad_u) + tau_p = 0")
problem.add_equation("dt(u)  + grad(p) + cross(fz_hat, u) - nu*div(grad_u) + lift(tau2u, -1) = -dot(u0,grad(u)) - dot(u,grad(u0))")

problem.add_equation("integ(p) = 0")
if (isNoSlip):
    # no-slip BCs
    problem.add_equation("dot(u, ex)(x='left')  = 0")
    problem.add_equation("dot(u, ex)(x='right')  = 0")
    problem.add_equation("dot(u, ey)(x='left')  = 0")
    problem.add_equation("dot(u, ey)(x='right') = 0")
    problem.add_equation("dot(u, ez)(x='left')  = 0")
    problem.add_equation("dot(u, ez)(x='right') = 0")
else:
    # stress-free BCs
    problem.add_equation("dot(u, ex)(x='left')      = 0")
    problem.add_equation("dot(u, ex)(x='right')     = 0")
    problem.add_equation("dot(dx(u), ey)(x='left')  = 0")
    problem.add_equation("dot(dx(u), ey)(x='right') = 0")
    problem.add_equation("dot(dx(u), ez)(x='left')  = 0")
    problem.add_equation("dot(dx(u), ez)(x='right') = 0")

dom = u.domain
local_slice = dom.dist.grid_layout.slices(dom,scales=1)
gshape = dom.dist.grid_layout.global_shape(dom,scales=1)
u_array = np.zeros((N, 3) + np.shape(np.zeros(gshape)[local_slice]))

cp_filenames = glob(load_cp)
# cp_indices = [int(fname.split('checkpoint_s')[-1][:-3]) for fname in glob("/home3/loconno2/mhd/ecs5/data/1/checkpoint/*")]
solver = problem.build_solver(d3.RK222)


for cp_index in range(len(cp_filenames)):
    with h5py.File(cp_filenames[cp_index - 1], "r") as file:
        # logger.info(cp_filenames[cp_index - 1])
        uarr_index = round(N * file['scales']['sim_time'][()][0] / T)        
        if uarr_index >= N:
            continue
        # u_array[uarr_index, ...] = file['tasks']['u'][()][0, ...][local_slice]
        # logger.info(np.shape(file['tasks']['u'][()][0, ...][(:,) + local_slice]))
        uy_local = file['tasks']['u'][()][0, ...][0, ...][local_slice]
        uz_local = file['tasks']['u'][()][0, ...][1, ...][local_slice]
        ux_local = file['tasks']['u'][()][0, ...][2, ...][local_slice]
        # logger.info(np.shape(uy))
        u_array[uarr_index, 0, ...] = uy_local.copy()
        u_array[uarr_index, 1, ...] = uz_local.copy()
        u_array[uarr_index, 2, ...] = ux_local.copy()

third = np.prod(gshape)
vecSize = np.prod(gshape)*3
vec = np.ones(vecSize)
def vecToField(solver,vec):
    solver.state[1]['g']
    solver.state[1].change_scales(1)
    solver.state[1]['g'][0] = vec[:third].reshape(gshape)[local_slice]
    solver.state[1]['g'][1] = vec[third:2*third].reshape(gshape)[local_slice]
    solver.state[1]['g'][2] = vec[2*third:3*third].reshape(gshape)[local_slice]

def fieldToVec(solver):
    vecu = np.zeros(gshape)
    vecv = np.zeros(gshape)
    vecw = np.zeros(gshape)
    solver.state[1]['g']
    solver.state[1].change_scales(1)
    vecu[local_slice] = solver.state[1]['g'][0]
    vecv[local_slice] = solver.state[1]['g'][1]
    vecw[local_slice] = solver.state[1]['g'][2]
    vecu = CW.allreduce(vecu,op=MPI.SUM).reshape(third)
    vecv = CW.allreduce(vecv,op=MPI.SUM).reshape(third)
    vecw = CW.allreduce(vecw,op=MPI.SUM).reshape(third)
    vec = np.hstack((vecu,vecv,vecw))
    return vec

def monodromyMult(q0):
    global solver, sim_index, timestep, T, N
    # problem.time['g'] = 0 # Reset time
    # solver.stop_sim_time = T - timestep
    start_time = solver.sim_time

    vecToField(solver,q0)
    try:
        logger.info('simulation index = {}'.format(sim_index))
        sim_index += 1
        for step_index in range(N):
            u0.change_scales(1)
            u0['g'] = u_array[step_index, ...]
            solver.step(timestep)
            # if (solver.iteration-1) % 100 == 0:
            #     logger.info('Iteration=%i, Time=%e, dt=%e' %(solver.iteration, solver.sim_time, timestep))
            # step_index += 1
    except:
        logger.error('Exception raised, triggering end of main loop.')
        raise
    qT = fieldToVec(solver)
    end_time = solver.sim_time
    # logger.info(end_time - start_time)
    return qT

import scipy.sparse as sp
from scipy import linalg
logger.info('Running first simulation to construct Linear Operator')
Psi = sp.linalg.LinearOperator((vecSize,vecSize),matvec=monodromyMult)
logger.info('Running additional simulations to satisfy EVP tolerance...')
mu, v = sp.linalg.eigs(Psi,k=NEV,tol=evp_tol)
CW.barrier()
if(CW.rank==0):
    print(mu)
    # np.savez(fileName,eigs=mu,modes=v)
    # np.savez(fileName+'_eigsOnly',eigs=mu)
