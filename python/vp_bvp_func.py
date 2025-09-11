import numpy as np
import dedalus.public as d3
from mpi4py import MPI
CW = MPI.COMM_WORLD
import logging
logger = logging.getLogger(__name__)
import sys

def vp_bvp_func(binit):
    domain = binit.domain
    dist = domain.dist
    bases = domain.bases
    Fbases = [basis for basis in bases if isinstance(basis, d3.RealFourier)]
    xbasis = [basis for basis in bases if not isinstance(basis, d3.RealFourier)][0]
    # print(dist.coords)
    # sys.exit()
    # Bases
    coords = d3.CartesianCoordinates('y', 'z', 'x')
    # coords = domain.coords

    # Fields
    phi = dist.Field(name='phi', bases=bases)
    A = dist.VectorField(coords, name='A', bases=bases)
    b = dist.VectorField(coords, name='b', bases=bases)
    binit.change_scales(1)
    b['g'] = binit['g'].copy()

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
    bvp_problem = d3.LBVP(variables=[A, phi, tau1A, tauphi], namespace=locals())

    bvp_problem.add_equation((d3.trace(grad_A), 0))
    bvp_problem.add_equation((d3.curl(A) + grad_phi + lift(tau1A), b))

    bvp_problem.add_equation("Ay(x='left') = 0", condition="(ny!=0) or (nz!=0)")
    bvp_problem.add_equation("Az(x='left') = 0", condition="(ny!=0) or (nz!=0)")
    bvp_problem.add_equation("Ay(x='right') = 0", condition="(ny!=0) or (nz!=0)")
    bvp_problem.add_equation("Az(x='right') = 0", condition="(ny!=0) or (nz!=0)")

    bvp_problem.add_equation("Ax(x='left') = 0", condition="(ny==0) and (nz==0)")
    bvp_problem.add_equation("Ay(x='left') = 0", condition="(ny==0) and (nz==0)")
    bvp_problem.add_equation("Az(x='left') = 0", condition="(ny==0) and (nz==0)")
    bvp_problem.add_equation("phi(x='left') = 0", condition="(ny==0) and (nz==0)")

    # Build solver
    solver = bvp_problem.build_solver()
    solver.solve()
    # CW.barrier()
    # print(CW.rank)
    logger.info('bvp solved.')
    return A['g'].copy()
