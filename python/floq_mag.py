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
import csv

try:    
    args = docopt(__doc__)
    filename = Path(args['<config_file>'])
except:
    filename = path + "/mri_options.cfg"
config = ConfigParser()
config.read(str(filename))

logger.info(config.items('parameters'))

scale = config.getfloat('parameters','scale')
logger_cadence = config.getint('parameters','logger_cadence')
write_cadence = config.getint('parameters','write_cadence')
Nsims = config.getint('parameters','Nsims')

suffix = eval(config.get('parameters', 'suffix'))
load_cp = eval(config.get('parameters','load_cp'))
cp_filenames = glob(eval(config.get('parameters','cp_filenames')))

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
else:
    Ro = config.getfloat('parameters', 'Ro')
    Ra = config.getfloat('parameters', 'Ra')

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

timestep = T / N

phi = dist.Field(name='phi', bases=bases)
A = dist.VectorField(coords, name='A', bases=bases)
u0 = dist.VectorField(coords, name='u0', bases=bases)
fz_hat = dist.VectorField(coords, name='fz_hat', bases=xbasis)
fz_hat['g'][1] = f

t = dist.Field()
# tau_phi = dist.Field()
tauphi = dist.Field(name='tauphi')

tau1A = dist.VectorField(coords, name='tau1A', bases=(ybasis,zbasis))
tau2A = dist.VectorField(coords, name='tau2A', bases=(ybasis,zbasis))

lift_basis = xbasis.clone_with(a=1/2, b=1/2) # First derivative basis
lift = lambda A, n: d3.Lift(A, lift_basis, n)
grad_A = d3.grad(A) + ex*lift(tau1A,-1)

problem = d3.IVP([phi, A, tau1A, tau2A, tauphi], time=t, namespace=locals())
problem.add_equation("trace(grad_A) + tauphi = 0")
problem.add_equation("dt(A)  + grad(phi) - eta*div(grad_A) + lift(tau2A, -1) = 0")
# problem.add_equation("dt(A)  + grad(phi) - eta*div(grad_A) + lift(tau2A, -1) = cross(u0, curl(A))")
problem.add_equation("phi(x='left')  = 0")
problem.add_equation("phi(x='right') = 0")
problem.add_equation("integ(phi) = 0")

problem.add_equation("dot(A, ey)(x='left')  = 0")
problem.add_equation("dot(A, ez)(x='left')  = 0")
problem.add_equation("dot(A, ey)(x='right') = 0")
problem.add_equation("dot(A, ez)(x='right') = 0")

dom = A.domain
local_slice = dom.dist.grid_layout.slices(dom,scales=1)
gshape = dom.dist.grid_layout.global_shape(dom,scales=1)

scales = 1.0
local_slice_scaled = dom.dist.grid_layout.slices(dom,scales=scales)
gshape_scaled = dom.dist.grid_layout.global_shape(dom,scales=scales)
u_array = np.zeros((N, 3) + np.shape(np.zeros(gshape_scaled)[local_slice_scaled]))

linalg_scales = 1.0
local_slice_linalg = dom.dist.grid_layout.slices(dom,scales=linalg_scales)
gshape_linalg = dom.dist.grid_layout.global_shape(dom,scales=linalg_scales)

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
        uy_local = file['tasks']['u'][()][0, ...][0, ...][local_slice_scaled]
        uz_local = file['tasks']['u'][()][0, ...][1, ...][local_slice_scaled]
        ux_local = file['tasks']['u'][()][0, ...][2, ...][local_slice_scaled]
        # logger.info(np.shape(uy))
        u_array[uarr_index, 0, ...] = uy_local.copy()
        u_array[uarr_index, 1, ...] = uz_local.copy()
        u_array[uarr_index, 2, ...] = ux_local.copy()

third = np.prod(gshape_linalg)
vecSize = np.prod(gshape_linalg)*3
vec = np.ones(gshape_linalg)
def vecToField(solver,vec):
    solver.state[1]['g']
    solver.state[1].change_scales(linalg_scales)
    solver.state[1]['g'][0] = vec[:third].reshape(gshape_linalg)[local_slice_linalg]
    solver.state[1]['g'][1] = vec[third:2*third].reshape(gshape_linalg)[local_slice_linalg]
    solver.state[1]['g'][2] = vec[2*third:3*third].reshape(gshape_linalg)[local_slice_linalg]
    solver.state[1].change_scales(1)

def fieldToVec(solver):
    vecu = np.zeros(gshape_linalg)
    vecv = np.zeros(gshape_linalg)
    vecw = np.zeros(gshape_linalg)
    solver.state[1]['g']
    solver.state[1].change_scales(linalg_scales)
    vecu[local_slice_linalg] = solver.state[1]['g'][0]
    vecv[local_slice_linalg] = solver.state[1]['g'][1]
    vecw[local_slice_linalg] = solver.state[1]['g'][2]
    vecu = CW.allreduce(vecu,op=MPI.SUM).reshape(third)
    vecv = CW.allreduce(vecv,op=MPI.SUM).reshape(third)
    vecw = CW.allreduce(vecw,op=MPI.SUM).reshape(third)
    vec = np.hstack((vecu,vecv,vecw))
    return vec

def monodromyMult():
    global solver, timestep, T, N
    # solver.stop_sim_time = T - timestep
    # start_time = solver.sim_time
    # problem.time['g'] = 0 # Reset time
    # solver = problem.build_solver(d3.RK222)


    # vecToField(solver,q0)
    try:
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
    # qT = fieldToVec(solver)
    # end_time = solver.sim_time
    # logger.info(end_time - start_time)
    return
    # return qT

import scipy.sparse as sp
from scipy import linalg
# logger.info('Running first simulation to construct Linear Operator')
if load_cp == 'default':
    A.fill_random()
else:
    raise
    with open(load_cp, 'rb') as f:
        ic0 = np.load(f)
mag0 = ((A@A)**0.5).evaluate().allreduce_L2_norm(normalize_volume=True)
ev_list = []
for ii in range(Nsims):
    # logger.info(mag0)
    monodromyMult()
    mag0 = ((A@A)**0.5).evaluate().allreduce_L2_norm(normalize_volume=True)
    ev_list.append(np.log(mag0))
    logger.info('sim_index={}; ev={}'.format(ii, np.log(mag0)))
    if ii % write_cadence == 0 and CW.rank == 0:
        with open(path + '/data/ev_sequence.npy', 'wb') as f:
            # wr = csv.writer(f, quoting=csv.QUOTE_ALL)
            np.save(f, np.array(ev_list))
            # wr.writerow(ev_list)
        # with open(path + "/data/mode" + str(ii).zfill(4) + ".npy", 'wb') as f:
        #     np.save(f, ic0)

# Psi = sp.linalg.LinearOperator((vecSize,vecSize),matvec=monodromyMult)
# logger.info('Running additional simulations to satisfy EVP tolerance...')
# mu, v = sp.linalg.eigs(Psi,k=NEV,tol=evp_tol)
# CW.barrier()
# if(CW.rank==0):
#     print(mu)
#     writeName = path + '/' + suffix + 'floquet'
#     np.savez(writeName,eigs=mu,modes=v)
#     np.savez(writeName+'_eigsOnly',eigs=mu)
