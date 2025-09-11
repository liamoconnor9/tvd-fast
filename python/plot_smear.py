from glob import glob
from docopt import docopt
from configparser import ConfigParser
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
plt.ioff()
from dedalus.extras import plot_tools
import logging
import sys
logger = logging.getLogger(__name__)
import os
path = os.path.dirname(os.path.abspath(__file__))
from mpi4py import MPI
CW = MPI.COMM_WORLD

def plot_shmear(filename, start, count):
    global sin_1d_lst, cos_1d_lst, indices, sim_times, x_g
    # Plot settings
    title_func = lambda sim_time: 't = {:.1f}'.format(sim_time)
    tasks = ['cos_1d_uy', 'sin_1d_uy']
    with h5py.File(filename, mode='r') as file:
        # for n, task in enumerate(tasks):
        # x_key = [key for key in file['scales'].keys() if 'x_hash' == key[:6]]
        # x_g = file['scales'][x_key]
        dset_sin = file['tasks']['sin_1d_uy']
        dset_cos = file['tasks']['cos_1d_uy']
        for index in range(start, start+count):
            full_index = file['scales/write_number'][index]
            sim_time_n = file['scales/sim_time'][index]
            indices.append(int(full_index))
            sim_times.append(sim_time_n)
            sin_1d_lst.append(dset_sin[index, 0, 0])
            cos_1d_lst.append(dset_cos[index, 0, 0])

if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync
    global suffix, ar, ary, arz, last_index, isHydro, solver_name, sin_1d_lst, cos_1d_lst, indices, sim_times, x_g

    if len(sys.argv) > 1:
        suffix = sys.argv[1]
        if suffix[-1] == '/':
            suffix = suffix[:-1]
    else:
        suffix = 'dnsNy4_T1e4'

    filename = "{}/{}/options.cfg".format(path, suffix)
    config = ConfigParser()
    config.read(str(filename))

    Ly = eval(config.get('parameters','Ly'))
    Lz = eval(config.get('parameters','Lz'))
    Lx = eval(config.get('parameters','Lx'))
    isHydro = config.getboolean('parameters','isHydro')
    solver_name = eval(config.get('parameters','SOLVER'))
    sin_1d_lst = []
    cos_1d_lst = []
    indices = []
    sim_times = []
    if CW.rank == 0:
        print(solver_name)
        print(solver_name == 'ql_fly.py')

    ary = Ly / Lx
    arz = Lz / Lx 

    slicepoints = glob("{}/{}/slicepoints/*.h5".format(path, suffix))
    last_index = len(slicepoints)
    with h5py.File(slicepoints[0], mode='r') as file:
        # for n, task in enumerate(tasks):
        x_key = [key for key in file['scales'].keys() if 'x_hash' == key[:6]][0]
        x_g = file['scales'][x_key][()]

    # Create output directory if needed
    output_mid_path_cleany=pathlib.Path('{}/{}/cleany'.format(path, suffix))
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_mid_path_cleany.exists():
                output_mid_path_cleany.mkdir()

    output_mid_path_cleanx=pathlib.Path('{}/{}/cleanx'.format(path, suffix))
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_mid_path_cleanx.exists():
                output_mid_path_cleanx.mkdir()

    output_mid_path_cleanz=pathlib.Path('{}/{}/cleanz'.format(path, suffix))
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_mid_path_cleanz.exists():
                output_mid_path_cleanz.mkdir()

    post.visit_writes(slicepoints, plot_shmear)
    sin_data = np.array(sin_1d_lst).T
    cos_data = np.array(cos_1d_lst).T
    # print(np.shape(sin_data))
    # print(np.array(sim_times) / sim_times[1])
    complex_shmear = cos_data + 1j*sin_data
    Nt = 10000
    sim_times_uniform = np.linspace(sim_times[0], sim_times[-1], Nt)
    uniform_shmear = np.zeros((np.shape(complex_shmear)[0], Nt), dtype=np.complex128)
    for row in range(np.shape(complex_shmear)[0]):
        uniform_shmear[row, :] = np.interp(sim_times_uniform, np.array(sim_times), complex_shmear[row, :])
    print(np.shape(uniform_shmear))
    # for 
    X, Y = np.meshgrid(sim_times_uniform, x_g)

    fig, axs = plt.subplots(2, sharex=True)
    fig.suptitle(r'$\langle u_y \exp(iz) \rangle_{y,z}$')

    axs[0].set_title('real')
    axs[0].set_ylabel(r'$x$')
    axs[0].pcolormesh(X, Y, uniform_shmear.real)
    axs[1].set_title('imag')
    axs[1].set_ylabel(r'$x$')
    axs[1].set_xlabel(r'$t$')
    axs[1].pcolormesh(X, Y, uniform_shmear.imag)
    # plt.pcolormesh(cos_data, aspect=100)
    savefigname = '{}/{}/shmear.png'.format(path, suffix)
    plt.savefig(savefigname)
    print(savefigname)
    plt.close()

    U, S, Vh = np.linalg.svd(uniform_shmear, full_matrices=False)
    print(np.allclose(uniform_shmear, np.dot(U * S, Vh)))

    plt.scatter(list(range(len(S))), S)
    plt.title('SVD values')
    savefigname = '{}/{}/svd_vals.png'.format(path, suffix)
    plt.savefig(savefigname)
    print(savefigname)

    plt.close()
    # rep = np.dot(U * S, Vh)
    
    S0 = np.zeros_like(S)
    S0[0] = S[0]
    mode0 = np.dot(U * S0, Vh)
    plt.close()

    fig, axs = plt.subplots(2, sharex=True)
    fig.suptitle('SVD mode 0')

    axs[0].set_title('real')
    axs[0].set_ylabel(r'$x$')
    axs[0].pcolormesh(X, Y, mode0.real)
    axs[1].set_title('imag')
    axs[1].set_ylabel(r'$x$')
    axs[1].set_xlabel(r'$t$')
    axs[1].pcolormesh(X, Y, mode0.imag)
    # plt.pcolormesh(cos_data, aspect=100)
    savefigname = '{}/{}/mode0.png'.format(path, suffix)
    plt.savefig(savefigname)
    print(savefigname)
    
    S1 = np.zeros_like(S)
    S1[1] = S[1]
    mode1 = np.dot(U * S1, Vh)
    plt.close()
    
    complex_tvec0 = Vh[0, :]
    plt.plot(sim_times_uniform, complex_tvec0.real, label='real')
    plt.plot(sim_times_uniform, complex_tvec0.imag, label='imag')
    plt.plot(sim_times_uniform, np.abs(complex_tvec0), label='mag')
    plt.legend()
    savefigname = '{}/{}/tvec0.png'.format(path, suffix)
    plt.savefig(savefigname)
    print(savefigname)

    plt.close()
    fig, axs = plt.subplots(2, sharex=True)
    fig.suptitle('SVD mode 1')

    axs[0].set_title('real')
    axs[0].set_ylabel(r'$x$')
    axs[0].pcolormesh(X, Y, mode1.real)
    axs[1].set_title('imag')
    axs[1].set_ylabel(r'$x$')
    axs[1].set_xlabel(r'$t$')
    axs[1].pcolormesh(X, Y, mode1.imag)
    # plt.pcolormesh(cos_data, aspect=100)
    savefigname = '{}/{}/mode1.png'.format(path, suffix)
    plt.savefig(savefigname)
    print(savefigname)
    plt.close()
    
    complex_tvec1 = Vh[1, :]
    plt.plot(sim_times_uniform, complex_tvec1.real, label='real')
    plt.plot(sim_times_uniform, complex_tvec1.imag, label='imag')
    plt.plot(sim_times_uniform, np.abs(complex_tvec1), label='mag')
    plt.legend()
    savefigname = '{}/{}/tvec1.png'.format(path, suffix)
    plt.savefig(savefigname)
    print(savefigname)

    plt.close()
    
    S2 = np.zeros_like(S)
    S1[2] = S[2]
    mode2 = np.dot(U * S2, Vh)
    fig, axs = plt.subplots(2, sharex=True)
    fig.suptitle('SVD mode 2')

    axs[0].set_title('real')
    axs[0].set_ylabel(r'$x$')
    axs[0].pcolormesh(X, Y, mode2.real)
    axs[1].set_title('imag')
    axs[1].set_ylabel(r'$x$')
    axs[1].set_xlabel(r'$t$')
    axs[1].pcolormesh(X, Y, mode2.imag)
    # plt.pcolormesh(cos_data, aspect=100)
    savefigname = '{}/{}/mode2.png'.format(path, suffix)
    plt.savefig(savefigname)
    print(savefigname)
    plt.close()
    
    complex_tvec2 = Vh[2, :]
    plt.plot(sim_times_uniform, complex_tvec2.real, label='real')
    plt.plot(sim_times_uniform, complex_tvec2.imag, label='imag')
    plt.plot(sim_times_uniform, np.abs(complex_tvec2), label='mag')
    plt.legend()
    savefigname = '{}/{}/tvec2.png'.format(path, suffix)
    plt.savefig(savefigname)
    print(savefigname)
