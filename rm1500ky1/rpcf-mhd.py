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
init_timestep = config.getfloat('parameters', 'init_timestep')
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

if 'floquet' in suffix:
    B0 = b * np.exp(growth_rate * t)
    A0 = A * np.exp(growth_rate * t)
if doLoadVelocity:
    vars = []
else:
    vars = [p, u, taup, tau1u, tau2u]
if not isHydro:
    vars += [phi, A, tau1A, tau2A]
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


        if not isHydro and not isKinematic:
            if is2d:
                NS_RHS -= d3.Integrate((d3.cross(b, d3.curl(b))), 'y') / Ly
            else:
                NS_RHS -= d3.cross(b, d3.curl(b))
        problem.add_equation((NS_LHS,   NS_RHS))
        # boundary conditions
        problem.add_equation("integ(p)              = 0") 
        problem.add_equation("dot(u, ex)(x='left')  = 0")
        problem.add_equation("dot(u, ex)(x='right') = 0")
        problem.add_equation("dot(u, ey)(x='left')  = dot(U0, ey)(x='left')")
        problem.add_equation("dot(u, ey)(x='right') = dot(U0, ey)(x='right')")
        problem.add_equation("dot(u, ez)(x='left')  = dot(U0, ez)(x='left')")
        problem.add_equation("dot(u, ez)(x='right') = dot(U0, ez)(x='right')")

    # induction
    if not isHydro:
        IND_LHS = dt(A) + d3.grad(phi) - 1 / Rm * d3.div(grad_A) + lift(tau2A, -1)
        IND_RHS = d3.cross(u, b)
        problem.add_equation((IND_LHS, IND_RHS))
        DIVA_LHS = d3.trace(grad_A)
        DIVA_RHS = 0
        problem.add_equation((DIVA_LHS, DIVA_RHS))
        problem.add_equation("phi(x='left')  = 0")
        problem.add_equation("phi(x='right') = 0")
        problem.add_equation("dot(A, ey)(x='left')  = 0")
        problem.add_equation("dot(A, ez)(x='left')  = 0")
        problem.add_equation("dot(A, ey)(x='right') = 0")
        problem.add_equation("dot(A, ez)(x='right') = 0")



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
A0 = A.copy()
A1 = A.copy()
u0 = u.copy()
u1 = u.copy()

def inner_p(field1, field2):
    return integ(field1 @ field2)

def get_norm(field):
    return np.sqrt(inner_p(field, field))

def analyze_floquet_func():
    return

if load_cp == 'default':
    u.fill_random(seed=seed)
    u.low_pass_filter(scales=0.0625)
    u.change_scales(1)
    u['g'] *= 1e-5 * x * (Lx - x)
    u.change_scales(1)
    u['g'][0] += S * x
    if not isHydro:
        logger.info("populating magnetic potential with noisy Bz={} initial condition".format(B0_z))
        logger.info('populating velocity with noise initial condition')
        A.fill_random(seed=2*seed)
        A.low_pass_filter(scales=0.0625)
        A['g'] *= ic_scale_A

        # curlA = d3.Curl(A).evaluate()
        # curlA.change_scales(1)
        # curlA['g'] *= 1e-6 * (x - Lx/2) * (x + Lx/2)
        # A.change_scales(1)
        # A['g'] = vp_bvp_func(curlA.copy())

        A['g'][0] += -2*B0_z*np.cos(np.pi*x / Lx) / (np.pi / Lx)
        A['g'] *= (x - Lx/2) * (x + Lx/2)
elif load_cp == 'noise':
    u.fill_random(seed=seed)
    u.low_pass_filter(scales=0.0625)
    u.change_scales(1)
    u['g'] *= 1e-5 * x * (Lx - x)
    u.change_scales(1)
    u['g'][0] += S * x
    if not isHydro:
        logger.info("populating magnetic potential with noisy Bz={} initial condition".format(B0_z))
        logger.info('populating velocity with noise initial condition')
        A.fill_random(seed=seed*2)
        # A.low_pass_filter(scales=0.0625)
        A['g'] *= 1e-5 * (x - Lx/2) * (x + Lx/2)
        
        
        # curlA = d3.Curl(A).evaluate()
        # curlA.change_scales(1)
        # curlA['g'] *= ic_scale_A * (x - Lx/2) * (x + Lx/2)
        # A.change_scales(1)
        # A['g'] = vp_bvp_func(curlA.copy())
else:
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
        if not isHydro:
            try:
                A.load_from_hdf5(file, 0, task='A')
                A.change_scales(1)
                A['g'] *= ic_scale_A
            except:
                logger.info('failed to load vector potential (magnetic field) data. Continuing with just the flow state assuming we loaded from hydro...')
                # A['g'][0] = -2*B0_z*np.cos(np.pi*x / Lx) / (np.pi / Lx)
                A.fill_random()
                A.low_pass_filter(scales=0.0625)
                A['g'] *= ic_scale_A * (x - Lx/2) * (x + Lx/2)
                logger.info('appending noisy magnetic field to existing hydro initial condition')
        imported_time = file['scales']['sim_time'][()][0]
        if doLoadTimestep:
            init_timestep = file['scales']['timestep'][()][0]
    if 'kin' in suffix:
        logger.info('normalizing magnetic potential...')
        A_normed = (A / get_norm(A)).evaluate()
        A_normed.change_scales(1)
        A.change_scales(1)
        A['g'] = A_normed['g'].copy()

# u.change_scales(1)
A.change_scales(1)
if doLoadVelocity:
    uinit.change_scales(dealias)
uinit['g'] = u['g'].copy()
Ainit['g'] = A['g'].copy()
binit.change_scales(dealias)
binit['g'] = b.evaluate()['g'].copy()

useCFL = True
if not doLoadVelocity and not isKinematic:
    CFL = d3.CFL(solver, initial_dt=init_timestep, cadence=10, safety=cfl_safety, threshold=0.05,
             max_change=1.5, min_change=0.5, max_dt=max_timestep)
    CFL.add_velocity(u)
    logger.info('constructing CFL with velocity')
    if not isKinematic and not isHydro:
        CFL.add_velocity(b)
        logger.info('appending magnetic field to CFL criteria')

elif not isKinematic and not isHydro:
    logger.info('constructing CFL for use with only magnetic field')
    CFL = d3.CFL(solver, initial_dt=init_timestep, cadence=10, safety=cfl_safety, threshold=0.05,
             max_change=1.5, min_change=0.5, max_dt=max_timestep)
    CFL.add_velocity(b)
else:
    useCFL = False
    timestep = init_timestep

mode_sin_lst = []
mode_cos_lst = []
Nmodes_track = 10
for i in range(Nmodes_track):
    ki = i + 1
    mode_temp = dist.Field(name='mode_sin{}'.format(ki), bases=zbasis)
    mode_temp['g'] = np.sin(ki * 2 * np.pi * z / Lz)
    mode_sin_lst.append(mode_temp.copy())
    mode_temp = dist.Field(name='mode_cos{}'.format(ki), bases=zbasis)
    mode_temp['g'] = np.cos(ki * 2 * np.pi * z / Lz)
    mode_cos_lst.append(mode_temp.copy())

def get_z_mode_amplitude(mode_number, field):
    mode_sin = mode_sin_lst[mode_number - 1]
    mode_cos = mode_cos_lst[mode_number - 1]
    projected_field = np.sqrt(integz(field * mode_sin)**2) + np.sqrt(integz(field * mode_cos)**2)
    return integy(integx((projected_field))) / vol


ux_1d = integy(integx(u @ ex)) / Ly / Lx
bdotgrad_bx_1d = integy(integx((b @ d3.grad(b)) @ ex)) / Ly / Lx

A_record = A
b_record = b
u_record = u
if 'floquet' in suffix:
    A_record = A0
    b_record = B0
    u_record = u * y_unity
elif is2d:
    u_record = u * y_unity


sparse_output = isKinematic and not isHydro
# sparse_output = isKinematic

if scalars_sim_dt != 0:
    scalars = solver.evaluator.add_file_handler(path + '/scalars', sim_dt=scalars_sim_dt, max_writes=1000, mode=fh_mode)
    scalars.add_task(d3.Integrate(0.5 * (u_record @ u_record)) / vol, name = 'ke')
    scalars.add_task(d3.Integrate(0.5 * (u_record @ ey)**2) / vol, name = 'ke_y')
    scalars.add_task(d3.Integrate(0.5 * (u_record @ ez)**2) / vol, name = 'ke_z')
    scalars.add_task(d3.Integrate(0.5 * (u_record @ ex)**2) / vol, name = 'ke_x')
    scalars.add_task(d3.Integrate((u_record - uinit) @ (u_record - uinit)) / vol, name = 'udiff')
    scalars.add_task(d3.Integrate(u_record @ d3.curl(u_record) / vol), name = 'uhelicity')
    if not isHydro:
        num_keff = d3.Integrate(d3.curl(b_record) @ d3.curl(b_record))
        den_keff = d3.Integrate(b_record @ b_record)
        scalars.add_task(np.sqrt(num_keff / den_keff), name = 'keff')

        scalars.add_task(d3.Integrate(A_record @ b_record / vol), name = 'bhelicity')
        scalars.add_task(d3.Integrate((A_record - Ainit) @ (A_record - Ainit)) / vol, name = 'Adiff')
        scalars.add_task(d3.Integrate(0.5 * (b_record @ b_record)) / vol, name = 'be')
        scalars.add_task(d3.Integrate(0.5 * (b_record @ ey)**2) / vol, name = 'be_y')
        scalars.add_task(d3.Integrate(0.5 * (b_record @ ez)**2) / vol, name = 'be_z')
        scalars.add_task(d3.Integrate(0.5 * (b_record @ ex)**2) / vol, name = 'be_x')
        if not sparse_output:
            B0_profileyz = integy(integz(b_record)) / Ly /Lz
            B0_pertyz = b_record - B0_profileyz
            B0_profileyz_energy = integx(0.5 * B0_profileyz @ B0_profileyz) / Lx
            B0_pertyz_energy = integ(0.5 * B0_pertyz @ B0_pertyz) / Ly / Lz / Lx
            scalars.add_task(B0_profileyz_energy, name='B0_profileyz_energy')
            scalars.add_task(B0_pertyz_energy, name='B0_pertyz_energy')
            scalars.add_task(B0_pertyz_energy / B0_profileyz_energy, name='B0_energy_ratio_yz')

            B0_profiley = integy(b_record) / Ly
            B0_perty = b_record - B0_profiley
            B0_profiley_energy = integz(integx(0.5 * B0_profiley @ B0_profiley)) / Lz / Lx
            B0_perty_energy = integ(0.5 * B0_perty @ B0_perty)  / Ly / Lz / Lx
            scalars.add_task(B0_profiley_energy, name='B0_profiley_energy')
            scalars.add_task(B0_perty_energy, name='B0_perty_energy')
            scalars.add_task(B0_perty_energy / B0_profiley_energy, name='B0_energy_ratio_y')

            norm_Ainit = get_norm(Ainit)
            Ainit_normed = Ainit / norm_Ainit
            norm_A = get_norm(A)
            A_normed = A / norm_A

            norm_binit = get_norm(binit)
            binit_normed = binit / norm_binit
            norm_b = get_norm(b)
            b_normed = b / norm_b

            projA = inner_p(Ainit_normed, A_normed)
            projb = inner_p(binit_normed, b_normed)
            scalars.add_task(projA, name='projA')
            scalars.add_task(projb, name='projb')

    def add_mode_n(n):
        scalars.add_task(get_z_mode_amplitude(n, u@u), name = 'ke_mode{}'.format(n))
        if not isHydro:
            scalars.add_task(get_z_mode_amplitude(n, b_record@b_record), name = 'be_mode{}'.format(n))
    if not sparse_output:
        for ki in range(Nmodes_track):
            add_mode_n(ki + 1)

mode2_sin = np.sin(2*z_field)
mode2_cos = np.cos(2*z_field)

mode1_sinz = np.sin(z_field)
mode1_cosz = np.cos(z_field)

mode1_siny = np.sin(y_field)
mode1_cosy = np.cos(y_field)

if sp_sim_dt != 0 :
    slicepoints = solver.evaluator.add_file_handler(path + '/slicepoints', sim_dt=sp_sim_dt, max_writes=50, mode=fh_mode)

    if not sparse_output:
        slicepoints.add_task(integz(integy((u @ ey)*mode2_sin)), name = "sin_1d_uy")
        slicepoints.add_task(integz(integy((u @ ey)*mode2_cos)), name = "cos_1d_uy")

        slicepoints.add_task(integz(integy((b @ ey)*mode2_sin)), name = "sin_1d_by")
        slicepoints.add_task(integz(integy((b @ ey)*mode2_cos)), name = "cos_1d_by")

        slicepoints.add_task(integz(integy((b @ ey)*mode1_siny*mode1_sinz)), name = "sinysinz_by")
        slicepoints.add_task(integz(integy((b @ ey)*mode1_siny*mode1_cosz)), name = "sinycosz_by")
        slicepoints.add_task(integz(integy((b @ ey)*mode1_cosy*mode1_cosz)), name = "cosycosz_by")
        slicepoints.add_task(integz(integy((b @ ey)*mode1_cosy*mode1_sinz)), name = "cosysinz_by")

        slicepoints.add_task(integz(integy((u @ ey)*mode1_siny*mode1_sinz)), name = "sinysinz_uy")
        slicepoints.add_task(integz(integy((u @ ey)*mode1_siny*mode1_cosz)), name = "sinycosz_uy")
        slicepoints.add_task(integz(integy((u @ ey)*mode1_cosy*mode1_cosz)), name = "cosycosz_uy")
        slicepoints.add_task(integz(integy((u @ ey)*mode1_cosy*mode1_sinz)), name = "cosysinz_uy")

        slicepoints.add_task(u_record @ d3.curl(u_record)(y = 'center'), name = "{}_mid{}".format('uhelicity_midy', 'y'), scales=sp_scale)
        slicepoints.add_task(u_record @ d3.curl(u_record)(z = 'center'), name = "{}_mid{}".format('uhelicity_midz', 'z'), scales=sp_scale)
        slicepoints.add_task(u_record @ d3.curl(u_record)(x = 'center'), name = "{}_mid{}".format('uhelicity_midx', 'x'), scales=sp_scale)
        # slicepoints.add_task(u @ d3.curl(u)(x = 'center'), name = "{}_mid{}".format('uhelicity_midx', 'x'))

    vector_slices_tuples = [(u_record, 'v'), (d3.curl(u_record), 'omega')]
    if not isHydro:
        vector_slices_tuples += [(b_record, 'b'), (b_record @ d3.grad(b_record), 'b.grad_b')]
        if not sparse_output:
            slicepoints.add_task((A_record @ b_record)(y = 'center'), name = 'bhelicity_midy', scales=sp_scale)
            slicepoints.add_task((A_record @ b_record)(z = 'center'), name = 'bhelicity_midz', scales=sp_scale)
            slicepoints.add_task((A_record @ b_record)(x = 'center'), name = 'bhelicity_midx', scales=sp_scale)

    for field, field_name in vector_slices_tuples:
        for d2, unit_vec in zip(('x', 'y', 'z'), (ex, ey, ez)):
            slicepoints.add_task(d3.dot(field, unit_vec)(x = 'center'), name = "{}{}_mid{}".format(field_name, d2, 'x'), scales=sp_scale)
            slicepoints.add_task(d3.dot(field, unit_vec)(y = 'center'), name = "{}{}_mid{}".format(field_name, d2, 'y'), scales=sp_scale)
            slicepoints.add_task(d3.dot(field, unit_vec)(z = 'center'), name = "{}{}_mid{}".format(field_name, d2, 'z'), scales=sp_scale)

            # slicepoints.add_task(d3.Integrate(d3.dot(field, unit_vec), 'x') / Lx, name = "{}{}_avg{}".format(field_name, d2, 'x'), scales=sp_scale)
            # slicepoints.add_task(d3.Integrate(d3.dot(field, unit_vec), 'y') / Ly, name = "{}{}_avg{}".format(field_name, d2, 'y'), scales=sp_scale)
            # slicepoints.add_task(d3.Integrate(d3.dot(field, unit_vec), 'z') / Lz, name = "{}{}_avg{}".format(field_name, d2, 'z'), scales=sp_scale)

    # if not isHydro:
    #     for d2, unit_vec in zip(('x', 'y', 'z'), (ex, ey, ez)):
    #         if not is2d:
    #             slicepoints.add_task(d3.dot(d3.dot(b_record, d3.grad(b_record)), unit_vec)(y='center'), name='b.grad_b{}_midy'.format(d2))
    #         slicepoints.add_task(d3.dot(integy(d3.dot(b_record, d3.grad(b_record))) / Ly, unit_vec), name='b.grad_b{}_avgy'.format(d2))

if cp_sim_dt != 0:
    checkpoint = solver.evaluator.add_file_handler(path + '/checkpoint', max_writes=1, sim_dt=cp_sim_dt, mode=fh_mode)
    checkpoint.add_task(u, name = 'u', layout='g', scales=cp_scale)
    if not isHydro:
        checkpoint.add_task(A_record, name = 'A', layout='g', scales=cp_scale)
        checkpoint.add_task(b_record, name = 'b', layout='g', scales=cp_scale)


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
    # if 'kin' in suffix or 'floquet' in suffix:
    #     analyze_floquet_func()
    
    if useCFL:
        timestep = CFL.compute_timestep()
    if (solver.iteration-1) % logger_cadence == 0:
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
        if not isHydro:
            max_bdotb = flow.max('b.b')
            norm_b = np.sqrt(flow.volume_integral('b.b') / Ly / Lz / Lx)
            max_bdotgrad_b_meany = flow.max('b.grad_b_meany')
            relative_bdotgrad_b_meany = max_bdotgrad_b_meany / max_bdotb
            mean_A_norm = np.sqrt(flow.volume_integral('A_norm'))
            max_Rm = flow.max('Rm')
            loop_message += "max(Rm)={}; ".format(max_Rm)
            loop_message += "A_norm={}; ".format(mean_A_norm)
            loop_message += "relative_bdotgrad_b_meany={}; ".format(relative_bdotgrad_b_meany)
            loop_message += "norm_b={}; ".format(norm_b)
        logger.info(loop_message)
    solver.step(timestep)