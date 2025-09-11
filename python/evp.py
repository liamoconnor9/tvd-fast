"""
3D cartesian MRI initial value problem using vector potential formulation
Usage:
    mri.py <config_file> <seed>
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
import matplotlib.pyplot as plt

try:    
    args = docopt(__doc__)
    filename = Path(args['<config_file>'])
    seed = int(args['<seed>'])
except:
    logger.warning('invalid config file supplied. using default... w/ seed = 1')
    filename = path + "/mri_options.cfg"
    seed = 1

config = ConfigParser()
config.read(str(filename))

logger.info(config.items('parameters'))

scale = config.getfloat('parameters','scale')
ra_scale = config.getfloat('parameters','ra_scale')

suffix = eval(config.get('parameters', 'suffix'))
logger.info(suffix)
load_cp = eval(config.get('parameters','load_cp'))
solve_EVP = config.getboolean('parameters','solve_EVP')
logger_cadence = config.getint('parameters','logger_cadence')

Ny = config.getint('parameters','Ny')
Ly = eval(config.get('parameters','Ly'))

Nz = config.getint('parameters','Nz')
Lz = eval(config.get('parameters','Lz'))

Nx = config.getint('parameters','Nx')
Lx = eval(config.get('parameters','Lx'))

ic_scale = config.getfloat('parameters','ic_scale')
tau = config.getfloat('parameters','tau')
nu = config.getfloat('parameters','nu')
Pm = config.getfloat('parameters','Pm')
eta = nu / Pm

Ro = config.getfloat('parameters', 'Ro')
f = config.getfloat('parameters', 'f')
S = -1
f = -S / Ro
# S = -f * Ro

B0_coeff = config.getfloat('parameters', 'B0_coeff')
B0_z = np.sqrt( -(S * f * Lx) / B0_coeff ) / np.pi
# B0_z = B0_coeff * (-f * Lx**2 * S * np.pi**(-2))

logger.info('[INPUTTED] Ro = {}'.format(Ro))
logger.info('[INPUTTED] B0_coeff = {}'.format(B0_coeff))

logger.info('[OUTPUTTED] S = {}'.format(S))
logger.info('[OUTPUTTED] B0_z = {}'.format(B0_z))

# sys.exit()

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

cp_sim_dt = config.getfloat('parameters','cp_sim_dt')
cp_scale = config.getfloat('parameters','cp_scale')
sp_sim_dt = config.getfloat('parameters','sp_sim_dt')

ncpu = MPI.COMM_WORLD.size
log2 = np.log2(ncpu)
if log2 == int(log2) and not solve_EVP:
    mesh = [int(2**np.ceil(log2/2)),int(2**np.floor(log2/2))]
elif solve_EVP:
    mesh = None
else:
    logger.error("pretty sure this shouldn't happen... log2(ncpu) is not an int?")
    
logger.info("running on processor mesh={}".format(mesh))

# Bases
coords = d3.CartesianCoordinates('y', 'z', 'x')
dealias = 3/2

Nky = 22
Nkz = 22

ky_vec_global = np.linspace(0.1, 1, Nky)
kz_vec_global = np.linspace(0.1, 1, Nkz)

ky_mat_global = np.zeros((Nky, Nkz))
kz_mat_global = np.zeros((Nky, Nkz))
for row in range(Nky):
    ky_mat_global[row, :] = ky_vec_global[row]

for col in range(Nkz):
    kz_mat_global[:, col] = kz_vec_global[col]

omega_mat_local = np.zeros_like(ky_mat_global, dtype=np.complex128)
for row in range(Nky):
    for col in range(Nkz):
        if Nkz*row + col % CW.size != CW.rank:
            continue
        ky = ky_mat_global[row, col]
        kz = kz_mat_global[row, col]

        dist = d3.Distributor(coords, dtype=np.complex128, comm=MPI.COMM_SELF)
        Ny = Nz = 2
        Ly = 2 * np.pi / ky
        Lz = 2 * np.pi / kz
        ybasis = d3.ComplexFourier(coords['y'], size=Ny, bounds=(0, Ly))
        zbasis = d3.ComplexFourier(coords['z'], size=Ny, bounds=(0, Lz))
        xbasis = d3.ChebyshevT(coords['x'], size=Nx, bounds=(-Lx / 2.0, Lx / 2.0))

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
        B0 = dist.VectorField(coords, name='B0', bases=xbasis)

        loaded = np.load("/home3/loconno2/mhd/taun4_pm25_nu1en3/data/1/write_dict.npy", allow_pickle=True)[()]
        B0['g'][1] = B0_z
        # B0['g'][1] = 0.01*np.sin(np.pi*x/Lx)
        # print(loaded.keys())
        # print(np.shape(loaded['bz_avg']))

        # B0['g'][1] = B0_z

        U0 = dist.VectorField(coords, name='U0', bases=xbasis)
        U0['g'][0] = S * x

        # B0 = 0
        # B0 = dist.VectorField(coords, name='B0', bases=xbasis)
        # B0['g'][1] = 0

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

        tau3u = dist.Field(name='tau3u', bases=(xbasis,))

        lift_basis = xbasis.clone_with(a=1/2, b=1/2) # First derivative basis
        lift = lambda A, n: d3.Lift(A, lift_basis, n)


        # operations
        # b.store_last = True
        # b = d3.Curl(A) + ex*lift(tau1A,-1)
        b = d3.Curl(A)
        integy = lambda A: d3.Integrate(A, 'y')
        integz = lambda A: d3.Integrate(A, 'z')
        integx = lambda A: d3.Integrate(A, 'x')
        integ = lambda A: integy(integz(integx(A)))

        dx = lambda A: d3.Differentiate(A, coords['x'])
        omega = dist.Field(name='omega')
        dt = lambda argy: -1j*omega*argy
        grad_u = d3.grad(u) + ex*lift(tau1u,-1) # First-order reduction
        grad_A = d3.grad(A) + ex*lift(tau1A,-1) # First-order reduction
        grad_b = d3.grad(b)
        # grad_B0 = d3.grad(B0)


        DIVU_LHS = d3.trace(grad_u) + taup
        DIVU_RHS = 0

        DIVA_LHS = d3.trace(grad_A)
        DIVA_RHS = 0

        NS_LHS = dt(u) - nu*d3.div(grad_u) + d3.grad(p) + d3.cross(fz_hat, u) + lift(tau2u,-1)
        NS_RHS = d3.cross(u, d3.curl(u)) - d3.cross(b, d3.curl(b))
        # NS_RHS = d3.dot(b, grad_b) - d3.dot(u, grad_u)

        NS_LHS += d3.cross(B0, d3.curl(d3.curl(A))) + d3.cross(d3.curl(A), d3.curl(B0))
        NS_LHS += d3.dot(u, d3.grad(U0)) + d3.dot(U0, grad_u)

        IND_LHS = dt(A) + d3.grad(phi) - eta*d3.div(grad_A) + lift(tau2A,-1) - d3.cross(U0, d3.curl(A)) - d3.cross(u, B0)
        IND_RHS = 0

        problem = d3.EVP([p, phi, u, A, taup, tau1u, tau2u, tau1A, tau2A], namespace=locals(), eigenvalue=omega)

        problem.add_equation((DIVU_LHS, 0))
        problem.add_equation((DIVA_LHS, 0))
        problem.add_equation((NS_LHS,   0))
        problem.add_equation((IND_LHS,  0))

        if (isNoSlip):
            # no-slip BCs
            problem.add_equation("u(x='left')  = 0")
            problem.add_equation("u(x='right') = 0")
        else:
            # stress-free BCs
            problem.add_equation("dot(u, ex)(x='left')      = 0")
            problem.add_equation("dot(u, ex)(x='right')     = 0")
            problem.add_equation("dot(dx(u), ey)(x='left')  = 0")
            problem.add_equation("dot(dx(u), ey)(x='right') = 0")
            problem.add_equation("dot(dx(u), ez)(x='left')  = 0")
            problem.add_equation("dot(dx(u), ez)(x='right') = 0")

        # problem.add_equation("integ(integ(u,'x'),'y')@ey = 0") 
        problem.add_equation("integ(p)       = 0") 
        problem.add_equation("phi(x='left')  = 0")
        problem.add_equation("phi(x='right') = 0")

        problem.add_equation("dot(A, ey)(x='left')  = 0")
        problem.add_equation("dot(A, ez)(x='left')  = 0")
        problem.add_equation("dot(A, ey)(x='right') = 0")
        problem.add_equation("dot(A, ez)(x='right') = 0")


        solver = problem.build_solver(entry_cutoff=0)
        # for i,sp in enumerate(solver.subproblems):
        #     print(i)
        #     print(sp.group)

        solver.solve_dense(solver.subproblems[3])
        # solver.solve_sparse(solver.subproblems[3], 10, target=0.0)
        evals = solver.eigenvalues[np.isfinite(solver.eigenvalues)]
        evals = evals[np.argsort(evals.imag)]
        # print(evals[0])
        # plt.scatter(evals.real, evals.imag)
        # print(evals.imag)
        # sys.exit()
        # figname = "{}/c_plane_oishi.png".format(path)
        # plt.savefig(figname)
        # print(figname)
        # omega_mat_local[row, col] = ky*kz*1j
        omega_mat_local[row, col] = evals[-1]


CW.barrier()
if CW.rank == 0:
    CW.Reduce(MPI.IN_PLACE, omega_mat_local, op=MPI.SUM, root=0)
else:
    CW.Reduce(omega_mat_local, omega_mat_local, op=MPI.SUM, root=0)

CW.barrier()
if CW.rank == 0:
    plt.xlabel(r'$k_y$')
    plt.ylabel(r'$k_z$')
    plt.title('Growth Rate')
    plt.pcolor(ky_mat_global, kz_mat_global, omega_mat_local.imag)
    plt.colorbar()
    filename ="{}/spectrum.png".format(path)
    plt.savefig(filename)
    print(filename)
