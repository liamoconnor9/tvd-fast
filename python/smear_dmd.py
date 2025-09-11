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
from pydmd.plotter import plot_modes_2D
from pydmd import DMD, BOPDMD
from pydmd.plotter import plot_eigs, plot_summary
from pydmd.preprocessing import hankel_preprocessing

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
            if not full_index % 10 == 0:
                continue
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
    Nt = 1000
    sim_times_uniform = np.linspace(sim_times[0], sim_times[-1], Nt)
    uniform_shmear = np.zeros((np.shape(complex_shmear)[0], Nt), dtype=np.complex128)
    for row in range(np.shape(complex_shmear)[0]):
        uniform_shmear[row, :] = np.interp(sim_times_uniform, np.array(sim_times), complex_shmear[row, :])
    print(np.shape(uniform_shmear))
    # for 


    X, Y = np.meshgrid(sim_times_uniform, x_g)
    d = 3
    optdmd = BOPDMD(svd_rank=3, num_trials=2)

    # Wrap the model with the preprocessing routine.
    delay_optdmd = hankel_preprocessing(optdmd, d=d)

    # Fit the model to the noisy data.
    # Note: BOPDMD models need the data X and the times of data collection t for fitting.
    # Hence if we apply time-delay, we must adjust the length of our time vector accordingly.
    delay_t = sim_times_uniform[: -d + 1]
    delay_optdmd.fit(uniform_shmear.real, t=delay_t)

    # Plot a summary of the DMD results.
    plot_summary(delay_optdmd, x=x_g, d=d)

    # Print computed eigenvalues (frequencies are given by imaginary components).
    # Also plot the resulting data reconstruction.
    print(
        f"Frequencies (imaginary component): {np.round(delay_optdmd.eigs, decimals=3)}"
    )
    plt.title("Reconstructed Data")
    plt.imshow(delay_optdmd.reconstructed_data.real)
    savefigname = '{}/{}/dmd0.png'.format(path, suffix)
    plt.savefig(savefigname)
    print(savefigname)

