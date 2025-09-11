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

try:    
    args = docopt(__doc__)
    filename = Path(args['<config_file>'])
except:
    filename = path + "/options.cfg"

config = ConfigParser()
config.read(str(filename))
logger.info(config.items('parameters'))

suffix = eval(config.get('parameters', 'suffix'))
logger_cadence = config.getint('parameters','logger_cadence')
try:
    seed = config.getint('parameters','seed')
except:
    seed = 1
load_cp = eval(config.get('parameters','load_cp'))

try:
    is2d = config.getboolean('parameters', 'is2d')
except:
    is2d = False
isHydro = config.getboolean('parameters','isHydro')
isKinematic = config.getboolean('parameters','isKinematic')
try:
    doLoadVelocity = config.getboolean('parameters','doLoadVelocity')
except:
    doLoadVelocity = False
if doLoadVelocity and not isKinematic:
    raise

doLoadTimestep = config.getboolean('parameters','doLoadTimestep')


Ny = config.getint('parameters','Ny')
Ly = eval(config.get('parameters','Ly'))
Nz = config.getint('parameters','Nz')
Lz = eval(config.get('parameters','Lz'))
Nx = config.getint('parameters','Nx')
Lx = eval(config.get('parameters','Lx'))
# if is2d:
#     Ny = 4

Npx = 16
Npz = 64
# marginz = Lz / Npz
# marginx = Lx / Npx
delz = Lz / (Npz + 1)
delx = Lx / (Npx + 1)
zlocations = [(i + 1)*delz          for i in range(Npz)]
xlocations = [(i + 1)*delx - Lx / 2 for i in range(Npx)]
locations = []
for col in range(Npx):
    for row in range(Npz):
        locations.append((xlocations[col], zlocations[row]))
Np = Npx * Npz

# logger.info(xlocations[0])
# logger.info(xlocations[-1])
# logger.info(zlocations[0])
# logger.info(zlocations[-1])
# sys.exit()

ic_scale_u = config.getfloat('parameters','ic_scale_u')
ic_scale_A = config.getfloat('parameters','ic_scale_A')
Ro = config.getfloat('parameters','Ro')
Re = config.getfloat('parameters','Re')
Rm = config.getfloat('parameters','Rm')
B0_coeff = config.getfloat('parameters', 'B0_coeff')
try:
    growth_rate = config.getfloat('parameters', 'growth_rate')
except:
    growth_rate = -0.01952041
timestep = init_timestep = config.getfloat('parameters', 'init_timestep')
max_timestep = config.getfloat('parameters', 'max_timestep')
cfl_safety = config.getfloat('parameters', 'cfl_safety')
stop_sim_time = config.getfloat('parameters', 'stop_sim_time') + max_timestep
# wall_time = 60. * 60. * config.getfloat('parameters', 'wall_time_hr')
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
if S == 0:
    B0_z = B0_coeff
logger.info('[INPUTTED] Ro = {}'.format(Ro))
logger.info('[INPUTTED] B0_coeff = {}'.format(B0_coeff))
logger.info('[OUTPUTTED] S = {}'.format(S))
logger.info('[OUTPUTTED] B0_z = {}'.format(B0_z))

ncpu = MPI.COMM_WORLD.size
log2 = np.log2(ncpu)
mesh = [1, round(ncpu / 1)]

# log2 = np.log2(ncpu)
# if log2 == int(log2):
#     mesh = [int(2**np.ceil(log2/2)),int(2**np.floor(log2/2))]
# else:
#     logger.info("pretty sure this shouldn't happen... log2(ncpu) is not an int?")
#     mesh=None
    
logger.info("running on processor mesh={}".format(mesh))

# Bases
coords = d3.CartesianCoordinates('y', 'z', 'x')
dealias = 3/2
# dist = d3.Distributor(coords, dtype=np.float64)
dist = d3.Distributor(coords, dtype=np.float64, mesh=mesh)
ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(0, Ly), dealias=dealias)
zbasis = d3.RealFourier(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)
xbasis = d3.ChebyshevT(coords['x'], size=Nx, bounds=(-Lx / 2.0, Lx / 2.0), dealias=dealias)
y = dist.local_grid(ybasis)
z = dist.local_grid(zbasis)
x = dist.local_grid(xbasis)
ey = dist.VectorField(coords, name='ey')
ez = dist.VectorField(coords, name='ez')
ex = dist.VectorField(coords, name='ex')
ey['g'][0] = 1
ez['g'][1] = 1
ex['g'][2] = 1

bases = (ybasis,zbasis,xbasis)
tau_bases = (ybasis,zbasis)
if is2d:
    twod_bases = (zbasis,xbasis)
    twod_tau_bases = (zbasis, )
else:
    twod_bases = bases
    twod_tau_bases = tau_bases

if is2d:
    def integy(A):
        try:
            return d3.Integrate(A, 'y')
        except:
            return Ly*A
    
    def dy(A):
        try:
            return d3.Differentiate(A, coords['y'])
        except:
            return 0*A
    
else:
    integy = lambda A: d3.Integrate(A, 'y')
    dy = lambda A: d3.Differentiate(A, coords['y'])


integz = lambda A: d3.Integrate(A, 'z')
integx = lambda A: d3.Integrate(A, 'x')

integ = lambda A: integy(integz(integx(A)))
dz = lambda A: d3.Differentiate(A, coords['z'])
dx = lambda A: d3.Differentiate(A, coords['x'])
from dedalus.core.operators import TimeDerivative
dt = lambda argy: TimeDerivative(argy)
lift_basis = xbasis.clone_with(a=1/2, b=1/2) # First derivative basis
lift = lambda A, n: d3.Lift(A, lift_basis, n)

#List of vectors, fields, vector fields indexed by particle
particlelocations = []
particlevelocities = []
for ii in range(Np):
    particlelocations.append(dist.VectorField(coords,name='pl'+str(ii)))
    particlevelocities.append(dist.VectorField(coords,name='pv'+str(ii)))

def particle_gutter():
    if CW.rank == 0:
        for ii in range(Np):
            pre_x = particlelocations[ii]['g'][0, ...][0][0][0]
            if pre_x < 0:
                particlelocations[ii]['g'][0, ...][0][0][0] += Lx
                logger.info('PARTICLE OUT OF DOMAIN')
                # raise
            if pre_x > Lx:
                particlelocations[ii]['g'][0, ...][0][0][0] -= Lx
                logger.info('PARTICLE OUT OF DOMAIN')
                # raise

# nccs
U0 = dist.VectorField(coords, name='U0', bases=xbasis)
U0['g'][0] = S * x
z_hat = dist.VectorField(coords, name='z_hat', bases=xbasis)
z_hat['g'][1] = 1

t = dist.Field(name='t')
y_unity = dist.Field(name='y_unity', bases=(ybasis,))
y_field = dist.Field(name='y_field', bases=(ybasis,))
y_field['g'] = y
y_unity['g'] = 1.0
z_field = dist.Field(name='z_field', bases=(zbasis,))
z_field['g'] = z
p = dist.Field(name='p', bases=twod_bases)
u = dist.VectorField(coords, name='u', bases=twod_bases)
if doLoadVelocity:
    logger.info('locking velocity to grid b/c it is prescribed with doLoadVelocity=True')
    u = d3.Grid(u).evaluate()
uinit = dist.VectorField(coords, name='uinit', bases=twod_bases)
taup = dist.Field(name='taup')
tau1u = dist.VectorField(coords, name='tau1u', bases=twod_tau_bases)
tau2u = dist.VectorField(coords, name='tau2u', bases=twod_tau_bases)
grad_u = d3.grad(u) + ex*lift(tau1u,-1)

A = dist.VectorField(coords, name='A', bases=bases)
Ainit = dist.VectorField(coords, name='Ainit', bases=bases)
binit = dist.VectorField(coords, name='binit', bases=bases)
phi = dist.Field(name='phi', bases=bases)
tau1A = dist.VectorField(coords, name='tau1A', bases=tau_bases)
tau2A = dist.VectorField(coords, name='tau2A', bases=tau_bases)
grad_A = d3.grad(A) + ex*lift(tau1A,-1)
b = d3.Curl(A)

def update_velocities():

    for ii in range(Np):
        if CW.rank == 0:
            pre_z = particlelocations[ii]['g'][1, ...][0][0][0]
            pre_x = particlelocations[ii]['g'][2, ...][0][0][0]
            
        else:
            pre_z = 0
            pre_x = 0
        pre_z = CW.bcast(pre_z, root=0)
        pre_x = CW.bcast(pre_x, root=0)

        interp = u(x=pre_x)(z=pre_z).evaluate()
        particlevelocities[ii]['g'] = interp['g'].copy()

for ii in range(len(locations)):
    particlelocations[ii]['g'][0] = Ly / 2
    particlelocations[ii]['g'][1] = locations[ii][1]
    particlelocations[ii]['g'][2] = locations[ii][0]

vars = [p, u, taup, tau1u, tau2u] + particlelocations
try:
    problem = d3.IVP(vars, time=t, namespace=locals())

    # divergence-free velocity
    if not doLoadVelocity:
        DIVU_LHS = d3.trace(grad_u) + taup
        DIVU_RHS = 0
        problem.add_equation((DIVU_LHS, DIVU_RHS))
        # problem.add_equation("dt(amp) - 0.01952041*amp = 0")


        # incompressible momentum
        NS_LHS = dt(u) + d3.grad(p) + 1 / Ro * d3.cross(z_hat, u) - 1 / Re * d3.div(grad_u) + lift(tau2u,-1)
        NS_RHS = d3.cross(u, d3.curl(u))

        problem.add_equation((NS_LHS,   NS_RHS))
        # boundary conditions
        problem.add_equation("integ(p)              = 0") 
        problem.add_equation("dot(u, ex)(x='left')  = 0")
        problem.add_equation("dot(u, ex)(x='right') = 0")
        problem.add_equation("dot(u, ey)(x='left')  = dot(U0, ey)(x='left')")
        problem.add_equation("dot(u, ey)(x='right') = dot(U0, ey)(x='right')")
        problem.add_equation("dot(u, ez)(x='left')  = dot(U0, ez)(x='left')")
        problem.add_equation("dot(u, ez)(x='right') = dot(U0, ez)(x='right')")

    for ii in range(Np):
        problem.add_equation((dt(particlelocations[ii]), particlevelocities[ii])) 

    solver = problem.build_solver(timestepper)
except Exception as e:
    logger.info(e)
    sys.exit()
    
solver.stop_sim_time = stop_sim_time

# Initial conditions
pert_scales = 0.25
fh_mode = 'overwrite'
imported_time = 0.0

N_checkpoints = 0
u0 = u.copy()
u1 = u.copy()

def inner_p(field1, field2):
    return integ(field1 @ field2)

def get_norm(field):
    return np.sqrt(inner_p(field, field))

load_path = "{}/{}".format(path, load_cp).replace(suffix + "/", "")
if not load_path[-3:] == '.h5':
    cps_all = glob('{}/checkpoint/*h5'.format(load_path))
    indices = [int(cp_path.split('checkpoint_s')[-1][:-3]) for cp_path in cps_all]
    cps_sorted = [x for _, x in sorted(zip(indices, cps_all))]
    load_path = cps_sorted[-1]
    logger.info('grabbing last checkpoint...')
logger.info('loading checkpoint: {}'.format(load_path))

# solver.load_state(load_path)    
with h5py.File(load_path, "r") as file:
    u.load_from_hdf5(file, 0, task='u')
    u.change_scales(1)
    u['g'] *= ic_scale_u
    imported_time = file['scales']['sim_time'][()][0]
    if doLoadTimestep:
        init_timestep = file['scales']['timestep'][()][0]

uinit['g'] = u['g'].copy()
Ainit['g'] = A['g'].copy()
binit.change_scales(dealias)
binit['g'] = b.evaluate()['g'].copy()

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=logger_cadence)

flow.add_property(d3.dot(u,u)*Re, name='Re')
flow.add_property(0.5*d3.dot(u,u), name='Ke')
if not isHydro:
    flow.add_property(d3.dot(u,u)*Rm, name='Rm')
    flow.add_property((d3.dot(A,A)), name='A_norm')
    flow.add_property(b@b, name='b.b')
    flow.add_property(integy(b @ d3.grad(b)) / Ly, name='b.grad_b_meany')

# flow.add_property(d3.dot(b, ez)**2/2, name='Be_z')
flow.add_property(0.5*(d3.curl(u)@ey)**2, name='enstr_y')
flow.add_property(0.5*(d3.curl(u)@ez)**2, name='enstr_z')
flow.add_property(0.5*(d3.curl(u)@ex)**2, name='enstr_x')
flow.add_property(u@u, name='u.u')
solver.evaluator.evaluate_handlers((flow.properties, ))

logger.info('Starting main loop')

while solver.proceed:
    update_velocities()
    # particle_gutter()

    if True:
    # if (solver.iteration-1) % logger_cadence == 0:
        max_Re = flow.max('Re')
        mean_Ke = flow.grid_average('Ke')
        mean_enstr_y = flow.grid_average('enstr_y')
        mean_enstr_z = flow.grid_average('enstr_z')
        mean_enstr_x = flow.grid_average('enstr_x')
        mean_unorm = np.sqrt(flow.volume_integral('u.u') / Ly / Lz / Lx)

        stop = False
        stop = stop or np.isnan(max_Re)
        stop = stop or np.isnan(mean_Ke)
        stop = stop or np.isnan(mean_Ke)
        if not isHydro:
            stop = stop or np.isnan(flow.grid_average('A_norm'))

        stop = stop or timestep < 1e-6
        if stop:
            logger.info('something is NAN. Terminating simulation. Get your shit together.')
            sys.exit()
        loop_message = ""
        loop_message += "Iteration={}; ".format(solver.iteration)
        loop_message += "Time={}; ".format(solver.sim_time)
        loop_message += "dt={}; ".format(timestep)
        loop_message += "mean(unorm)={}; ".format(mean_unorm)
        loop_message += "max(Re)={}; ".format(max_Re)
        loop_message += "avg(Ke)={}; ".format(mean_Ke)
        flowline = ""
        if CW.rank == 0:
            for ball_index in range(3):
                locations_temp = particlelocations[ball_index]['g'].transpose()[0][0][0]
                if np.any(np.isnan(locations_temp)):
                    raise
                else:
                    loop_message += ", position{}={}".format(ball_index, str(locations_temp))

            #         flowline += "{},{},".format(locations_temp[0], locations_temp[1])
            #         # flowline += str(solver.sim_time)
            # with open('flow_toys.csv', 'a') as flowfile:
            #     flowfile.write(flowline[:-1] + '\n')

        logger.info(loop_message)
    solver.step(timestep)