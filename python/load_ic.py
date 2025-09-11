from configparser import ConfigParser
import numpy as np
import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
import h5py
import dedalus.public as d3
from mpi4py import MPI

CW = MPI.COMM_WORLD
import logging
logger = logging.getLogger(__name__)

def load_ic(u, b, filename, unitvectors, seed=42):

    config = ConfigParser()
    config.read(str(filename))

    scale = config.getfloat('parameters','scale')

    Ny = config.getint('parameters','Ny')
    ny = round(scale * Ny)
    Ly = eval(config.get('parameters','Ly'))

    Nz = config.getint('parameters','Nz')
    nz = round(scale * Nz)
    Lz = eval(config.get('parameters','Lz'))

    Nx = config.getint('parameters','Nx')
    nx = round(scale * Nx)
    Lx = eval(config.get('parameters','Lx'))

    by0 = config.getfloat('parameters','by0')
    bz0 = config.getfloat('parameters','bz0')
    bx0 = config.getfloat('parameters','bx0')

    uy0 = config.getfloat('parameters','uy0')
    uz0 = config.getfloat('parameters','uz0')
    ux0 = config.getfloat('parameters','ux0')

    dtype = np.float64
    ncpu = MPI.COMM_WORLD.size
    log2 = np.log2(ncpu)
    if log2 == int(log2):
        mesh = [int(2**np.ceil(log2/2)),int(2**np.floor(log2/2))]
    else:
        logger.error("pretty sure this shouldn't happen... log2(ncpu) is not an int?")

    # Bases
    coords = d3.CartesianCoordinates('y', 'z', 'x')
    dist = d3.Distributor(coords, dtype=dtype, mesh=mesh)
    dealias = 3/2

    xbasis = d3.ChebyshevT(coords['x'], size=nx, bounds=(-Lx / 2.0, Lx / 2.0), dealias=dealias)
    ybasis = d3.RealFourier(coords['y'], size=ny, bounds=(0, Ly), dealias=dealias)
    zbasis = d3.RealFourier(coords['z'], size=nz, bounds=(0, Lz), dealias=dealias)

    # Fields
    U = dist.VectorField(coords, name='u', bases=(ybasis,zbasis,xbasis))
    B = dist.VectorField(coords, name='b', bases=(ybasis,zbasis,xbasis))

    U.fill_random(layout='c', distribution='uniform', seed=seed)
    B.fill_random(layout='c', distribution='uniform', seed=int(seed + 1e6))

    u.change_scales(scale)
    b.change_scales(scale)

    u['g'] = U['g'].copy()
    b['g'] = B['g'].copy()

    u.change_scales(1)
    b.change_scales(1)

    vol = Ly * Lz * Lx
    ey, ez, ex = unitvectors

    integ = lambda A: d3.Integrate(d3.Integrate(d3.Integrate(A, 'y'), 'z'), 'x')
    
    u.fill_random(layout='c', seed=seed)
    b.fill_random(layout='c', seed=int(seed + 1e6))
    u.low_pass_filter(scales=scale)
    b.low_pass_filter(scales=scale)

    evaledu = integ(u).evaluate()['g']
    evaledb = integ(b).evaluate()['g']

    if (CW.rank == 0):
        umean = np.squeeze(evaledu / vol).ravel()
        bmean = np.squeeze(evaledb / vol).ravel()
    else:
        umean = np.zeros((3, ))
        bmean = np.zeros((3, ))

    CW.Bcast([umean, MPI.DOUBLE], root=0)
    CW.Bcast([bmean, MPI.DOUBLE], root=0)

    # u['g'] -= umean[:, np.newaxis, np.newaxis, np.newaxis]
    # b['g'] -= bmean[:, np.newaxis, np.newaxis, np.newaxis]

    uog_eval = d3.Integrate(0.5*(u @ u)).evaluate()
    bog_eval = d3.Integrate(0.5*(b @ b)).evaluate()
    
    uyog_eval = d3.Integrate(0.5*(u @ ey)**2).evaluate()
    uzog_eval = d3.Integrate(0.5*(u @ ez)**2).evaluate()
    uxog_eval = d3.Integrate(0.5*(u @ ex)**2).evaluate()
    byog_eval = d3.Integrate(0.5*(b @ ey)**2).evaluate()
    bzog_eval = d3.Integrate(0.5*(b @ ez)**2).evaluate()
    bxog_eval = d3.Integrate(0.5*(b @ ex)**2).evaluate()

    if (CW.rank == 0):
        uog = uog_eval['g'][0] / vol
        bog = bog_eval['g'][0] / vol

        uyog = uyog_eval['g'][0] / vol
        uzog = uzog_eval['g'][0] / vol
        uxog = uxog_eval['g'][0] / vol
        byog = byog_eval['g'][0] / vol
        bzog = bzog_eval['g'][0] / vol
        bxog = bxog_eval['g'][0] / vol
    else:
        uog = 0.0
        bog = 0.0

        uyog = 0.0
        uzog = 0.0
        uxog = 0.0
        byog = 0.0
        bzog = 0.0
        bxog = 0.0
    
    uog = CW.bcast(uog, root=0)
    bog = CW.bcast(bog, root=0)

    u['c'] /= (uog)**0.5
    b['c'] /= (bog)**0.5

    u['c'] *= 0.0
    # b['c'] *= 1.0e4

    # uyog = CW.bcast(uyog, root=0)
    # uzog = CW.bcast(uzog, root=0)
    # uxog = CW.bcast(uxog, root=0)
    # byog = CW.bcast(byog, root=0)
    # bzog = CW.bcast(bzog, root=0)
    # bxog = CW.bcast(bxog, root=0)

    # u['g'][0] *= (uy0 / uyog)**0.5
    # u['g'][1] *= (uz0 / uzog)**0.5
    # u['g'][2] *= (ux0 / uxog)**0.5

    # b['g'][0] *= (by0 / byog)**0.5
    # b['g'][1] *= (bz0 / bzog)**0.5
    # b['g'][2] *= (bx0 / bxog)**0.5


    # from numpy.random import RandomState
    # # rstate = RandomState(seed)
    # print('hi')
    # a = []
    # for i in range(1000):
    #     a.append(np.random.rayleigh())
    # # print(rstate.rayleigh())
    # import matplotlib.pyplot as plt
    # plt.hist(a)
    # plt.savefig(path + '/hist.png')
    # logger.info(path + '/hist.png')

    # sys.exit()
    
    # bog_eval = d3.Integrate(0.5*(b @ b)).evaluate()

    # if (CW.rank == 0):
    #     bog = bog_eval['g'][0] / vol
    # else:
    #     bog = 0.0
    
    # bog = CW.bcast(bog, root=0)
    # logger.info(bog)



    return (u, b)



# else:
#     psi.fill_random(layout='c', seed=2*seed)
#     A.fill_random(layout='c', seed=int(2*seed + 1))
#     psi.low_pass_filter(scales=scale)
#     A.low_pass_filter(scales=scale)

#     evaledu = integ(d3.Curl(psi)).evaluate()['g']
#     evaledb = integ(d3.Curl(A)).evaluate()['g']

#     if (CW.rank == 0):
#         umean = np.squeeze(evaledu / vol).ravel()
#         bmean = np.squeeze(evaledb / vol).ravel()
#     else:
#         umean = np.zeros((3, ))
#         bmean = np.zeros((3, ))

#     CW.Bcast([umean, MPI.DOUBLE], root=0)
#     CW.Bcast([bmean, MPI.DOUBLE], root=0)

#     u['g'] -= umean[:, np.newaxis, np.newaxis, np.newaxis]
#     # b['g'] -= bmean[:, np.newaxis, np.newaxis, np.newaxis]
#     binit.change_scales(dealias)
#     binit['g'] = d3.Curl(A).evaluate()['g'].copy() - bmean[:, np.newaxis, np.newaxis, np.newaxis]
#     binit['g'][1] *= 4

#     logger.info("number of nonzero coefficients: {}".format(np.count_nonzero(binit['c'])))
#     logger.info("number of clearly nonzero coefficients: {}".format((binit['c'] > 1e-3).sum()))
#     logger.info("number of available coefficients: {}".format(np.size(binit['c'])))

#     binit.change_scales(1)
#     A.change_scales(1)
#     pre_func = binit['g'].copy()
#     A['g'] = vp_bvp_func(binit['g'].copy(), dist, bases, coords)
#     post_func_field = d3.Curl(A).evaluate()
#     post_func_field.change_scales(1)
#     post_func = post_func_field['g'].copy()

#     logger.info('vp_bvp_func successful: ' + str(np.allclose(pre_func, post_func)))

#     bog_eval = integ(0.5*(b @ b)).evaluate()
#     uog_eval = integ(0.5*(u @ u)).evaluate()

#     if (CW.rank == 0):
#         bog = bog_eval['g'][0] / vol
#         uog = uog_eval['g'][0] / vol
#     else:
#         bog = 0.0
#         uog = 0.0

#     bog = CW.bcast(bog, root=0)
#     uog = CW.bcast(uog, root=0)

#     # normalize magnetic energy, zero kinetic
#     u['g'] *= (uog)**(-0.5)
#     A['g'] *= (bog)**(-0.5)

#     from numpy.random import RandomState
#     rstate = RandomState(seed)
#     # be_mean = rstate.rayleigh(ra_scale)
#     be_mean = 3.0
#     ke_mean = 0.7

#     A['g'] *= (be_mean)
#     u['g'] *= (ke_mean)
#     logger.info('setting mean magnetic energy to: ' + str(be_mean))
#     logger.info('setting mean kinetic energy to: ' + str(ke_mean))

#     A.change_scales(1)
#     A['g'][0] *= 0.5**2 - x**2
#     A['g'][1] *= 0.5**2 - x**2
