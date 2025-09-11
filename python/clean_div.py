import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)


def clean_div(dist, bases, coords, u_data):

    ybasis = bases[0]
    zbasis = bases[1]
    xbasis = bases[2]

    y, z, x = dist.local_grids(ybasis, zbasis, xbasis)

    # Fields
    logger.info('divergence cleaning initiated...')
    ui = dist.VectorField(coords, name='ui', bases=bases)
    uc = dist.VectorField(coords, name='uc', bases=bases)
    ui['g'] = u_data.copy()
    pi = dist.Field(name='pi', bases=bases)
    tau_p = dist.Field(name='tau_p')

    # Problem
    problem = d3.LBVP([pi, uc, tau_p], namespace=locals())
    problem.add_equation("lap(pi) + tau_p = -div(ui)")
    problem.add_equation("uc - grad(pi) = ui")

    problem.add_equation("integ(pi) = 0")

    # Solver
    solver = problem.build_solver()
    solver.solve()
    logger.info('divergence cleaning successful')

    uc.change_scales(1)
    return uc['g'].copy()