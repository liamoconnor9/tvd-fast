from scipy import stats
from scipy.optimize import root, root_scalar, minimize_scalar
import numpy as np
import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
import h5py
import dedalus.public as d3
from mpi4py import MPI
CW = MPI.COMM_WORLD
rank = CW.Get_rank()
sys.path.append("..") # Adds higher directory to python modules path.
from vp_bvp_func import *
from glob import glob
import logging
logger = logging.getLogger(__name__)
from docopt import docopt
from pathlib import Path
from configparser import ConfigParser
import matplotlib.pyplot as plt
import gc

def run_sim(ky, Rm, config):
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
    Ly = 2*np.pi / ky
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
    # Rm = config.getfloat('parameters','Rm')
    # if Rm == None:
        # Rm = Rm_guess
    B0_coeff = config.getfloat('parameters', 'B0_coeff')
    try:
        growth_rate = config.getfloat('parameters', 'growth_rate')
    except:
        growth_rate = -0.01952041
    init_timestep = config.getfloat('parameters', 'init_timestep')
    timestep = max_timestep = config.getfloat('parameters', 'max_timestep')
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
    S = -1
    B0_z = B0_coeff * (-Ro * Lx**2 * S * np.pi**(-2))
    if S == 0:
        B0_z = B0_coeff

    ncpu = MPI.COMM_WORLD.size
    log2 = np.log2(ncpu)
    mesh = [1, round(ncpu / 1)]    

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
    fh_mode = 'overwrite'
    A0 = A.copy()

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
                # A['g'][0] = -2*B0_z*np.cos(np.pi*x / Lx) / (np.pi / Lx)
                A.fill_random()
                A.low_pass_filter(scales=0.0625)
                A['g'] *= ic_scale_A * (x - Lx/2) * (x + Lx/2)
        imported_time = file['scales']['sim_time'][()][0]
        if doLoadTimestep:
            init_timestep = file['scales']['timestep'][()][0]

    if 'kin' in suffix:
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
        if not isKinematic and not isHydro:
            CFL.add_velocity(b)

    elif not isKinematic and not isHydro:
        CFL = d3.CFL(solver, initial_dt=init_timestep, cadence=10, safety=cfl_safety, threshold=0.05,
                max_change=1.5, min_change=0.5, max_dt=max_timestep)
        CFL.add_velocity(b)
    else:
        useCFL = False
        timestep = init_timestep

    # Flow properties
    flow = d3.GlobalFlowProperty(solver, cadence=logger_cadence)

    flow.add_property(d3.dot(u,u)*Re, name='Re')
    flow.add_property(0.5*d3.dot(u,u), name='Ke')
    if not isHydro:
        flow.add_property(d3.dot(u,u)*Rm, name='Rm')
        flow.add_property((d3.dot(A,A)), name='A_norm')
        flow.add_property(0.5 * (b @ ey)**2, name='be_y')
        flow.add_property(b@b, name='b.b')
        flow.add_property(integy(b @ d3.grad(b)) / Ly, name='b.grad_b_meany')

    # flow.add_property(d3.dot(b, ez)**2/2, name='Be_z')
    flow.add_property(0.5*(d3.curl(u)@ey)**2, name='enstr_y')
    flow.add_property(0.5*(d3.curl(u)@ez)**2, name='enstr_z')
    flow.add_property(0.5*(d3.curl(u)@ex)**2, name='enstr_x')
    flow.add_property(u@u, name='u.u')
    solver.evaluator.evaluate_handlers((flow.properties, ))

    # logger.info('Starting main loop')

    be_y_lst = []
    time_lst = []
    be_setup = integ(0.5 * (b @ ey)**2) / Ly / Lz / Lx
    while solver.proceed:
        if solver.iteration % 100 == 0:
            CW.Barrier()
            be_delayed = be_setup.evaluate()
            CW.Barrier()
            if rank == 0:
                be_y_lst.append(be_delayed['g'][0][0][0])
            else:
                be_y_lst.append(0)
            
            time_lst.append(solver.sim_time)
            CW.Barrier()

        # if useCFL:
        #     timestep = CFL.compute_timestep()
        if logger_cadence != 0 and (solver.iteration-1) % logger_cadence == 0:
            # max_Re = flow.max('Re')
            # mean_Ke = flow.grid_average('Ke')
            # mean_enstr_y = flow.grid_average('enstr_y')
            # mean_enstr_z = flow.grid_average('enstr_z')
            # mean_enstr_x = flow.grid_average('enstr_x')
            # mean_unorm = np.sqrt(flow.volume_integral('u.u') / Ly / Lz / Lx)

            # stop = False
            # stop = stop or np.isnan(max_Re)
            # stop = stop or np.isnan(mean_Ke)
            # stop = stop or np.isnan(mean_Ke)
            # if not isHydro:
            #     stop = stop or np.isnan(flow.grid_average('A_norm'))

            # stop = stop or timestep < 1e-6
            # if stop:
            #     logger.info('something is NAN. Terminating simulation. Get your shit together.')
            #     sys.exit()
            loop_message = ""
            loop_message += "Iteration={}; ".format(solver.iteration)
            loop_message += "Time={}; ".format(solver.sim_time)
            # loop_message += "dt={}; ".format(timestep)
            # loop_message += "mean(unorm)={}; ".format(mean_unorm)
            # loop_message += "max(Re)={}; ".format(max_Re)
            # loop_message += "avg(Ke)={}; ".format(mean_Ke)
            # if not isHydro:
            #     max_bdotb = flow.max('b.b')
            #     norm_b = np.sqrt(flow.volume_integral('b.b') / Ly / Lz / Lx)
            #     max_bdotgrad_b_meany = flow.max('b.grad_b_meany')
            #     relative_bdotgrad_b_meany = max_bdotgrad_b_meany / max_bdotb
            #     mean_A_norm = np.sqrt(flow.volume_integral('A_norm'))
            #     max_Rm = flow.max('Rm')
            #     loop_message += "be_y={}; ".format(be_y_lst[-1])
            #     loop_message += "max(Rm)={}; ".format(max_Rm)
            #     loop_message += "A_norm={}; ".format(mean_A_norm)
            #     loop_message += "relative_bdotgrad_b_meany={}; ".format(relative_bdotgrad_b_meany)
            #     loop_message += "norm_b={}; ".format(norm_b)
            logger.info(loop_message)
        solver.step(timestep)
    
    Npts = len(time_lst)
    cutoff = int(Npts / 2)
    # del solver
    # del problem
    vars_local = locals()
    time_lst = time_lst[cutoff:]
    be_y_lst = be_y_lst[cutoff:]
    keys_to_del = []
    for key in vars_local.keys():
        if not (key == 'time_lst' or key == 'be_y_lst'):
            keys_to_del.append(key)
    for key in keys_to_del:
        del vars_local[key]
    del vars_local
    gc.collect()
    return (time_lst, be_y_lst)
