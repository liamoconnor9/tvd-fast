"""
3D cartesian MRI initial value problem using vector potential formulation
Usage:
    mri.py <config_file>
"""
from docopt import docopt
from configparser import ConfigParser
from pathlib import Path
import numpy as np
import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
import h5py
import dedalus.public as d3
from mpi4py import MPI
CW = MPI.COMM_WORLD
import logging
import pathlib
from glob import glob
logger = logging.getLogger(__name__)

def vp_bvp_func(bdata, dist, bases, coords):
    Fbases = [basis for basis in bases if isinstance(basis, d3.RealFourier)]
    xbasis = [basis for basis in bases if not isinstance(basis, d3.RealFourier)][0]

    # Fields
    phi = dist.Field(name='phi', bases=bases)
    A = dist.VectorField(coords, name='A', bases=bases)
    b = dist.VectorField(coords, name='b', bases=bases)
    b['g'] = bdata.copy()

    tauphi = dist.Field(name='tauphi', bases=Fbases)
    tau1A = dist.VectorField(coords, name='tau1A', bases=Fbases)

    ey = dist.VectorField(coords, name='ey')
    ez = dist.VectorField(coords, name='ez')
    ex = dist.VectorField(coords, name='ex')

    ey['g'][0] = 1
    ez['g'][1] = 1
    ex['g'][2] = 1

    Ay = A @ ey
    Az = A @ ez
    Ax = A @ ex

    lift_basis = xbasis.derivative_basis(1) # First derivative basis
    lift = lambda A: d3.Lift(A, lift_basis, -1)
    grad_A = d3.grad(A) + ex*lift(tau1A) # First-order reduction
    grad_phi = d3.grad(phi) + ex*lift(tauphi)

    # b = d3.Curl(A).evaluate()

    logger.info('solving bvp for vector potential A given b')
    problem = d3.LBVP(variables=[A, phi, tau1A, tauphi], namespace=locals())

    problem.add_equation((d3.trace(grad_A), 0))
    problem.add_equation((d3.curl(A) + grad_phi + lift(tau1A), b))

    problem.add_equation("Ay(x='left') = 0", condition="(ny!=0) or (nz!=0)")
    problem.add_equation("Az(x='left') = 0", condition="(ny!=0) or (nz!=0)")
    problem.add_equation("Ay(x='right') = 0", condition="(ny!=0) or (nz!=0)")
    problem.add_equation("Az(x='right') = 0", condition="(ny!=0) or (nz!=0)")

    problem.add_equation("Ax(x='left') = 0", condition="(ny==0) and (nz==0)")
    problem.add_equation("Ay(x='left') = 0", condition="(ny==0) and (nz==0)")
    problem.add_equation("Az(x='left') = 0", condition="(ny==0) and (nz==0)")
    problem.add_equation("phi(x='left') = 0", condition="(ny==0) and (nz==0)")

    # Build solver
    solver = problem.build_solver()
    solver.solve()
    logger.info('bvp solved.')
    A.change_scales(1)
    return A['g'].copy()

try:    
    args = docopt(__doc__)
    filename = Path(args['<config_file>'])
    seed = 1
    # seed = int(args['<seed>'])
except:
    logger.warning('invalid config file supplied. using default... w/ seed = 1')
    filename = path + "/mri_options.cfg"
    seed = 1

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

B0_coeff = config.getfloat('parameters', 'B0_coeff')
B0_z = B0_coeff * (-f * Lx**2 * S * np.pi**(-2))
if S == 0:
    B0_z = B0_coeff

logger.info('[INPUTTED] Ro = {}'.format(Ro))
logger.info('[INPUTTED] B0_coeff = {}'.format(B0_coeff))

logger.info('[OUTPUTTED] S = {}'.format(S))
logger.info('[OUTPUTTED] B0_z = {}'.format(B0_z))

isNoSlip = config.getboolean('parameters','isNoSlip')

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

# nccs
U0 = dist.VectorField(coords, name='U0', bases=xbasis)
U0['g'][0] = S * x

fz_hat = dist.VectorField(coords, name='fz_hat', bases=xbasis)
fz_hat['g'][1] = f

# damping timescale
Rayleigh = dist.Field(name='Rayleigh')

Rossby = dist.Field(name='Rossby')

if isContinuation:
    Rayleigh['g'] = Ra_vec[0]
    Rossby['g'] = Ro_vec[0]
else:
    Rayleigh['g'] = Ra
    Rossby['g'] = Ro

# Fields
t = dist.Field(name='t')
p = dist.Field(name='p', bases=bases)
u = dist.VectorField(coords, name='u', bases=bases)
uinit = dist.VectorField(coords, name='uinit', bases=bases)
Ainit = dist.VectorField(coords, name='Ainit', bases=bases)
psi = dist.VectorField(coords, name='psi', bases=bases)
taup = dist.Field(name='taup')
tau1u = dist.VectorField(coords, name='tau1u', bases=(ybasis,zbasis))
tau2u = dist.VectorField(coords, name='tau2u', bases=(ybasis,zbasis))

if not isHydro:
    binit = dist.VectorField(coords, name='binit', bases=bases)
    phi = dist.Field(name='phi', bases=bases)
    A = dist.VectorField(coords, name='A', bases=bases)
    tau1A = dist.VectorField(coords, name='tau1A', bases=(ybasis,zbasis))
    tau2A = dist.VectorField(coords, name='tau2A', bases=(ybasis,zbasis))
    b = d3.Curl(A)

if isConvecting:
    s = dist.Field(name='s', bases=bases)
    tau1s = dist.Field(name='tau1s', bases=(ybasis,zbasis))
    tau2s = dist.Field(name='tau2s', bases=(ybasis,zbasis))

lift_basis = xbasis.clone_with(a=1/2, b=1/2) # First derivative basis
lift = lambda A, n: d3.Lift(A, lift_basis, n)

# operations
# b.store_last = True
# b = d3.Curl(A) + ex*lift(tau1A,-1)
integy = lambda A: d3.Integrate(A, 'y')
integz = lambda A: d3.Integrate(A, 'z')
integx = lambda A: d3.Integrate(A, 'x')
integ = lambda A: integy(integz(integx(A)))

dx = lambda A: d3.Differentiate(A, coords['x'])
from dedalus.core.operators import TimeDerivative
dt = lambda argy: TimeDerivative(argy)
grad_u = d3.grad(u) + ex*lift(tau1u,-1)
if not isHydro:
    grad_A = d3.grad(A) + ex*lift(tau1A,-1)
    grad_b = d3.grad(b)
    DIVA_LHS = d3.trace(grad_A)
    DIVA_RHS = 0

if isConvecting:
    grad_s = d3.grad(s) + ex*lift(tau1s,-1)
# grad_B0 = d3.grad(B0)

DIVU_LHS = d3.trace(grad_u) + taup
DIVU_RHS = 0

b_dot_grad_b_coeff = int(not isKinematic)
logger.info("isKinematic = {}".format(isKinematic))
logger.info("Magnetic tension term coefficient = {}".format(b_dot_grad_b_coeff))

NS_LHS = dt(u) - nu*d3.div(grad_u) + d3.grad(p) + d3.cross(fz_hat, u) + lift(tau2u,-1)
if isConvecting:
    NS_LHS -= ex*s
if isHydro:
    NS_RHS = d3.cross(u, d3.curl(u))
else:
    NS_RHS = d3.cross(u, d3.curl(u)) - b_dot_grad_b_coeff * d3.cross(b, d3.curl(b))
    IND_LHS = dt(A) + d3.grad(phi) - eta*d3.div(grad_A) + lift(tau2A,-1) 

    IND_RHS = d3.cross(u, b)
    # IND_RHS = d3.cross(integy(u)/Ly, b)
    if tau != 0:
        omega = ((integ(A@A)/Ly/Lz/Lx) - 1) / tau
        IND_RHS -= omega * A

if not isConvecting:
    if isHydro:
        vars = [p, u, taup, tau1u, tau2u]
    else:
        vars = [p, phi, u, A, taup, tau1u, tau2u, tau1A, tau2A]
else:
    if isHydro:
        vars = [p, s, u, taup, tau1s, tau2s, tau1u, tau2u]
    else:
        vars = [p, phi, s, u, A, taup, tau1s, tau2s, tau1u, tau2u, tau1A, tau2A]

problem = d3.IVP(vars, time=t, namespace=locals())

problem.add_equation((DIVU_LHS, DIVU_RHS))
problem.add_equation((NS_LHS,   NS_RHS))
if not isHydro:
    problem.add_equation((DIVA_LHS, DIVA_RHS))
    problem.add_equation((IND_LHS,  IND_RHS))
    problem.add_equation("phi(x='left')  = 0")
    problem.add_equation("phi(x='right') = 0")
    problem.add_equation("dot(A, ey)(x='left')  = 0")
    problem.add_equation("dot(A, ez)(x='left')  = 0")
    problem.add_equation("dot(A, ey)(x='right') = 0")
    problem.add_equation("dot(A, ez)(x='right') = 0")

if (isNoSlip):
    # no-slip BCs
    problem.add_equation("dot(u, ex)(x='left')      = 0")
    problem.add_equation("dot(u, ex)(x='right')     = 0")
    problem.add_equation("dot(u, ey)(x='left')  = dot(U0, ey)(x='left')")
    problem.add_equation("dot(u, ey)(x='right') = dot(U0, ey)(x='right')")
    problem.add_equation("dot(u, ez)(x='left')  = dot(U0, ez)(x='left')")
    problem.add_equation("dot(u, ez)(x='right') = dot(U0, ez)(x='right')")
else:
    # stress-free BCs
    problem.add_equation("dot(u, ex)(x='left')      = dot(U0, ex)(x='left')")
    problem.add_equation("dot(u, ex)(x='right')     = dot(U0, ex)(x='right')")
    problem.add_equation("dot(dx(u), ey)(x='left')  = dot(dx(U0), ey)(x='left')")
    problem.add_equation("dot(dx(u), ey)(x='right') = dot(dx(U0), ey)(x='right')")
    problem.add_equation("dot(dx(u), ez)(x='left')  = dot(dx(U0), ez)(x='left')")
    problem.add_equation("dot(dx(u), ez)(x='right') = dot(dx(U0), ez)(x='right')")

# problem.add_equation("integ(integ(u,'x'),'y')@ey = 0") 
problem.add_equation("integ(p)       = 0") 

if isConvecting:
    problem.add_equation("dt(s) - Ra/Pr*(u@ex) - (nu/Pr) * div(grad_s) + lift(tau2s,-1)= -u@grad_s")
    problem.add_equation("s(x='right')=0")
    problem.add_equation("s(x='left')=0")

solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time


# Initial conditions
pert_scales = 0.25
fh_mode = 'overwrite'
imported_time = 0.0

# if load_cp == 'default':
#     logger.info('setting fresh new initial condition...')
#     logger.info('(w/ complementary fresh new initial condition smell)')

#     psi_pert_mag = 1e-6
#     A_pert_mag = 1e-6

#     # psi.change_scales(1)
#     # psi.fill_random()
#     # psi.low_pass_filter(scales=pert_scales)
#     # psi['g'][0] *= psi_pert_mag * np.cos(x*np.pi / Lx)
#     # psi['g'][1] *= psi_pert_mag * np.cos(x*np.pi / Lx)
#     # psi['g'][2] *= psi_pert_mag * np.cos(x*np.pi / Lx)
#     # temp = d3.Curl(psi).evaluate()
#     # temp.change_scales(1)
#     # u['g'] = temp['g'].copy()

#     if isHydro:
#         u.fill_random()
#         u['g'] *= 1e-3 * (x - Lx/2) * (x + Lx/2)
#     else:
#         u['g'] = 0

#     u.change_scales(1)
#     u['g'][0] += S * x
#     if not isHydro:
#         A.change_scales(1)
#         A.fill_random()
#         A.low_pass_filter(scales=pert_scales)
#         A['g'][0] *= A_pert_mag * np.cos(x*np.pi / Lx)
#         A['g'][1] *= A_pert_mag * np.cos(x*np.pi / Lx)
#         A['g'][2] *= A_pert_mag * np.cos(x*np.pi / Lx)
#         A['g'][0] += -2*B0_z*np.cos(np.pi*x / Lx) / (np.pi / Lx)

if load_cp == 'default':
    # u.fill_random()
    # u.low_pass_filter(scales=0.0625)
    # u['g'] *= 1e-5 * x * (Lx - x)
    u.change_scales(1)
    u['g'] = 0.0
    u['g'][0] += S * x
    if not isHydro:
        logger.info("populating magnetic potential with noisy Bz={} initial condition".format(B0_z))
        logger.info('populating velocity with noise initial condition')
        A.fill_random()
        A.low_pass_filter(scales=0.0625)
        A['g'] *= (x - Lx/2) * (x + Lx/2)
        curlA = d3.Curl(A).evaluate()
        curlA.change_scales(1)
        curlA['g'] *= 1e-6 * (x - Lx/2) * (x + Lx/2)
        A.change_scales(1)
        A['g'] = vp_bvp_func(curlA['g'].copy(), dist, bases, coords)
        A['g'][0] += -2*B0_z*np.cos(np.pi*x / Lx) / (np.pi / Lx)
        # A['g'] *= 1e-6 * (x - Lx/2) * (x + Lx/2)
        # A['g'] = curlA.evaluate()['g'].copy()
elif load_cp == 'noise':
    logger.info('populating velocity with noise initial condition')
    u.fill_random()
    u.low_pass_filter(scales=0.0625)
    u['g'] *= ic_scale_u * x * (Lx - x)
    u.change_scales(1)
    # u['g'] = 0.0
    u['g'][0] += S * x
    if not isHydro:
        logger.info('magnetic potential with noise initial condition')
        A.fill_random()
        A.low_pass_filter(scales=0.0625)
        A['g'] *= (x - Lx/2) * (x + Lx/2)
        curlA = d3.Curl(A).evaluate()
        curlA.change_scales(1)
        curlA['g'] *= ic_scale_A * (x - Lx/2) * (x + Lx/2)
        A.change_scales(1)
        A['g'] = vp_bvp_func(curlA['g'].copy(), dist, bases, coords)
        # A['g'] *= 1e-6 * (x - Lx/2) * (x + Lx/2)
        # A['g'] = curlA.evaluate()['g'].copy()

else:
    load_path = "{}/{}".format(path, load_cp).replace(suffix + "/", "")
    logger.info('loading checkpoint: {}'.format(load_path))

    # solver.load_state(load_path)    
    with h5py.File(load_path, "r") as file:
        # logger.info(np.shape(file['tasks']['u'][()][0, :, :, :171, :]))
        # sys.exit()
        # u.load_from_global_grid_data(file['tasks']['u'][()][0, :, :, :171, :])
        # u.fill_random()
        # u['g'] *= 1e-6
        # u_temp = u['g'].copy()
        u.load_from_global_grid_data(ic_scale_u * file['tasks']['u'][()][0, ...])
        u.change_scales(1)
        # u['g'] += u_temp
        if loadFromPert:
            new_u = (u + U0).evaluate()
            new_u.change_scales(1)
            u.change_scales(1)
            u['g'] = new_u['g'].copy()
        if not isHydro:
            try:
                A.load_from_global_grid_data(ic_scale_A * file['tasks']['A'][()][0, ...])
                A.change_scales(1)
                Ainit['g'] = A['g'].copy()
            except:
                logger.info('failed to load vector potential (magnetic field) data. Continuing with just the flow state assuming we loaded from hydro...')
                # A.change_scales(1)
                # A['g'][0] = -2*B0_z*np.cos(np.pi*x / Lx) / (np.pi / Lx)
                A.fill_random()
                A.low_pass_filter(scales=0.0625)
                A['g'] *= (x - Lx/2) * (x + Lx/2)
                curlA = d3.Curl(A).evaluate()
                curlA.change_scales(1)
                curlA['g'] *= ic_scale_A * (x - Lx/2) * (x + Lx/2)
                A.change_scales(1)
                A['g'] = vp_bvp_func(curlA['g'].copy(), dist, bases, coords)
                # logger.info('B0_z = {}'.format(B0_z))
                logger.info('appending noisy magnetic field to existing hydro initial condition')
        imported_time = file['scales']['sim_time'][()][0]
        if doLoadTimestep:
            init_timestep = file['scales']['timestep'][()][0]
        else:
            init_timestep = max_timestep

    if not isHydro and ic_scale_A == -1:
    
        A.change_scales(1)
        # A.fill_random()
        # A.low_pass_filter(scales=0.125)
        # A['g'][0] *= 1e-6 * np.cos(x*np.pi / Lx)
        # A['g'][1] *= 1e-6 * np.cos(x*np.pi / Lx)
        # A['g'][2] *= 1e-6 * np.cos(x*np.pi / Lx)

        A['g'][0] += -2*B0_z*np.cos(np.pi*x / Lx) / (np.pi / Lx)
        logger.info('B0_z = {}'.format(B0_z))
        logger.info('appending sinusoidal magnetic field to existing hydro initial condition')

    logger.info('successfully loaded data from checkpoint at sim_time = {}, timestep = {}'.format(imported_time, init_timestep))

logger.info('recording initial A and u data...')
if not  isHydro:
    Ainit['g'] = A['g'].copy()
uinit['g'] = u['g'].copy()


if (CW.rank == 0):
    if not os.path.exists(path + '/data/' + str(seed)):
        os.makedirs(path + '/data/' + str(seed))

if crumb_sim_dt != 0:
    crumb = solver.evaluator.add_file_handler(path + '/data/' + str(seed) + '/crumb', max_writes=1, sim_dt=crumb_sim_dt, mode=fh_mode, parallel="gather")
    crumb.add_task(u, name = 'u', layout='g', scales=crumb_scale)
    if not isHydro:
        crumb.add_task(A, name = 'A', layout='g', scales=crumb_scale)

if cp_sim_dt != 0:
    checkpoint = solver.evaluator.add_file_handler(path + '/data/' + str(seed) + '/checkpoint', max_writes=1, sim_dt=cp_sim_dt, mode=fh_mode)
    checkpoint.add_task(u, name = 'u', layout='g')
    if not isHydro:
        checkpoint.add_task(A, name = 'A', layout='g')

if scalars_sim_dt != 0:
    scalars = solver.evaluator.add_file_handler(path + '/data/' + str(seed) + '/scalars', sim_dt=scalars_sim_dt, max_writes=100, mode=fh_mode)


    scalars.add_task(d3.Integrate(0.5 * (u @ ey)**2) / vol, name = 'ke_y')
    scalars.add_task(d3.Integrate(0.5 * (u @ ez)**2) / vol, name = 'ke_z')
    scalars.add_task(d3.Integrate(0.5 * (u @ ex)**2) / vol, name = 'ke_x')
    scalars.add_task(d3.Integrate(d3.dot(u - uinit, u - uinit)/2), name = 'udiff')

    if not isHydro:
        scalars.add_task(d3.Integrate(0.5 * (b @ ey)**2) / vol, name = 'be_y')
        scalars.add_task(d3.Integrate(0.5 * (b @ ez)**2) / vol, name = 'be_z')
        scalars.add_task(d3.Integrate(0.5 * (b @ ex)**2) / vol, name = 'be_x')

        # measure the thing we want to measure
        scalars.add_task(d3.Integrate(0.5 * ((b @ ey)**2 - (b @ ez)**2)) / vol, name = 'be_y-be_z')
        scalars.add_task(d3.Integrate(d3.curl(b)) / vol, name = 'j')
        scalars.add_task(d3.Integrate(d3.dot(Ainit, A)), name = 'proj_A0')
        scalars.add_task(d3.Integrate(d3.dot(Ainit, A)) / (d3.Integrate(d3.dot(A, A))**0.5), name = 'proj_norm_A0')
        scalars.add_task(d3.Integrate(A@A)/Ly/Lz/Lx, name='AdotA_mean')
        scalars.add_task(d3.Integrate((A - Ainit)@(A - Ainit)/2), name = 'Adiff')
        if tau != 0:
            scalars.add_task(omega, name='omega')

    scalars.add_task(Rossby, name='Rossby')
    scalars.add_task(Rayleigh, name='Rayleigh')

if sp_sim_dt != 0 :
    try:
        # fh_mode = 'overwrite'
        slicepoints = solver.evaluator.add_file_handler(path + '/data/' + str(seed) + '/slicepoints', sim_dt=sp_sim_dt, max_writes=50, mode=fh_mode)
        if not isHydro:
            slices_tuples = [(b, 'b'), ((u), 'v'), (d3.curl(b), 'j')]
        else:
            slices_tuples = [((u), 'v')]

        for field, field_name in slices_tuples:
            for d2, unit_vec in zip(('x', 'y', 'z'), (ex, ey, ez)):
                slicepoints.add_task(d3.dot(field, unit_vec)(x = 'center'), name = "{}{}_mid{}".format(field_name, d2, 'x'))
                slicepoints.add_task(d3.dot(field, unit_vec)(y = 'center'), name = "{}{}_mid{}".format(field_name, d2, 'y'))
                slicepoints.add_task(d3.dot(field, unit_vec)(z = 'center'), name = "{}{}_mid{}".format(field_name, d2, 'z'))
                
                slicepoints.add_task(d3.Integrate(d3.dot(field, unit_vec), 'x'), name = "{}{}_avg{}".format(field_name, d2, 'x'))
                slicepoints.add_task(d3.Integrate(d3.dot(field, unit_vec), 'y'), name = "{}{}_avg{}".format(field_name, d2, 'y'))
                slicepoints.add_task(d3.Integrate(d3.dot(field, unit_vec), 'z'), name = "{}{}_avg{}".format(field_name, d2, 'z'))
                    
            slicepoints.add_task(d3.Integrate(d3.Integrate(d3.dot(field, ey), 'y'), 'z') / Ly / Lz, name = "{}{}_avg".format(field_name, 'y'))
            slicepoints.add_task(d3.Integrate(d3.Integrate(d3.dot(field, ez), 'y'), 'z') / Ly / Lz, name = "{}{}_avg".format(field_name, 'z'))


            

    except:
        logger.info("failed to add slicepoints file handler.")

CFL = d3.CFL(solver, initial_dt=init_timestep, cadence=10, safety=cfl_safety, threshold=0.05,
             max_change=1.5, min_change=0.5, max_dt=max_timestep)
CFL.add_velocity(u)
if not isKinematic and not isHydro:
    CFL.add_velocity(b)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=1)
flow.add_property(d3.dot(u,u)/nu, name='Re')
# flow.add_property(d3.dot(b,b)/nu, name='Rm')
flow.add_property(0.5*d3.dot(u,u), name='Ke')
flow.add_property((u - uinit)@(u - uinit)/2, name='udiff')

# flow.add_property(d3.dot(b, ez)**2/2, name='Be_z')
flow.add_property(0.5*(d3.curl(u)@ey)**2, name='enstr_y')
flow.add_property(0.5*(d3.curl(u)@ez)**2, name='enstr_z')
flow.add_property(0.5*(d3.curl(u)@ex)**2, name='enstr_x')
flow.add_property(u@u, name='u.u')
if not isHydro:
    flow.add_property((A - Ainit)@(A - Ainit)/2, name='Adiff')
    flow.add_property(0.5*d3.dot(b,b), name='Be')
    flow.add_property(A@A, name='A.A')

solver.evaluator.evaluate_handlers((flow.properties, ))

# print(solver.sim_time)
# sys.exit()

# Main loop
# try:
logger.info('Starting main loop')
Ro_t, Ra_t = Ro, Ra
while solver.proceed:
    timestep = CFL.compute_timestep()
    if isContinuation:
        param_index = np.argmax(t_vec > solver.sim_time)
        Ro_t = Ro_vec[param_index]
        Ra_t = int(isConvecting) * Ra_vec[param_index]
    if (solver.iteration-1) % logger_cadence == 0:
        max_Re = flow.max('Re')
        # max_Rm = flow.max('Rm')
        mean_Ke = flow.grid_average('Ke')
        # mean_Be_z = flow.grid_average('Be_z')
        mean_enstr_y = flow.grid_average('enstr_y')
        mean_enstr_z = flow.grid_average('enstr_z')
        mean_enstr_x = flow.grid_average('enstr_x')
        udiff = flow.volume_integral('udiff') / Ly / Lz / Lx
        stop = False

        if not isHydro:
            Adiff = flow.volume_integral('Adiff') / Ly / Lz / Lx
            mean_Be = flow.grid_average('Be')
            mean_AdotA = flow.volume_integral('A.A') / Ly / Lz / Lx
            mean_udotu = flow.volume_integral('A.A') / Ly / Lz / Lx
            stop = stop or np.isnan(mean_Be)
            # stop = stop or (mean_Be < 1e-7 and solver.sim_time - imported_time > 
            mean_Anorm = np.sqrt(flow.volume_integral('A.A') / Ly / Lz / Lx)

        mean_unorm = np.sqrt(flow.volume_integral('u.u') / Ly / Lz / Lx)
        # if (solver.iteration - 1) % (10 * logger_cadence) == 0:
        #     # u['g'] /= mean_unorm
        #     A['g'] /= mean_Anorm

        stop = stop or np.isnan(max_Re)
        stop = stop or np.isnan(mean_Ke)
        stop = stop or timestep < 1e-6
        if stop:
            logger.info('something is NAN. Terminating simulation. Get your shit together.')
            sys.exit()
        
        
        if tau != 0 and not isHydro:
            omeg_mag = omega.evaluate()
        if CW.rank == 0:
            # logger.info('Iteration=%i, Time=%e, dt=%e, max(Re)=%e, max(Rm)=%e, grid_average(Ke)=%e, grid_average(Be)=%e, grid_average(Be_z)=%e, Rossby=%e, Rayleigh=%e' %(solver.iteration, solver.sim_time, timestep, max_Re, max_Rm, mean_Ke, mean_Be, mean_Be_z, Ro_t, Ra_t))
            # logger.info('Iteration=%i, Time=%e, dt=%e, max(Re)=%e, max(Rm)=%e, grid_average(Ke)=%e, grid_average(Be)=%e, grid_average(Be_z)=%e, Rossby=%e, Rayleigh=%e' %(solver.iteration, solver.sim_time, timestep, max_Re, max_Rm, mean_Ke, mean_Be, mean_Be_z, Ro_t, Ra_t))
            loop_message = ""
            loop_message += "Iteration={}; ".format(solver.iteration)
            loop_message += "Time={}; ".format(solver.sim_time)
            loop_message += "dt={}; ".format(timestep)
            loop_message += "max(Re)={}; ".format(max_Re)
            loop_message += "avg(Ke)={}; ".format(mean_Ke)
            if not isHydro:
                loop_message += "avg(Be)={}; ".format(mean_Be)
            loop_message += "avg(enstr_y)={}; ".format(mean_enstr_y)
            loop_message += "avg(enstr_z)={}; ".format(mean_enstr_z)
            loop_message += "avg(enstr_x)={}; ".format(mean_enstr_x)
            if not isHydro:
                loop_message += "mean(unorm)={}; ".format(mean_unorm)
                loop_message += "mean(Anorm)={}; ".format(mean_Anorm)
                if tau != 0:
                    loop_message += "omega={}; ".format(omeg_mag['g'])
            try:
                loop_message += "Ro={}; ".format(Ro_t)
            except:
                pass
            try:
                loop_message += "Ra={}; ".format(Ra_t)
            except:
                pass
            try:
                loop_message += "udiff={}; ".format(udiff)
            except:
                pass
            try:
                loop_message += "Adiff={}; ".format(Adiff)
            except:
                pass
            logger.info(loop_message)
    U0.change_scales(1)
    U0['g'] = 0.0
    U0['g'][0] = -f * x
    U0['g'] *= Ro_t
    Rossby['g'] = Ro_t
    Rayleigh['g'] = Ra_t
    # CW.barrier()
    solver.step(timestep)