import numpy as np
# import dedalus.public as d3
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from mpi4py import MPI
import logging
logger = logging.getLogger(__name__)
CW = MPI.COMM_WORLD

def get_color(label):
    if label == 'be_y':
        return 'blue'
    elif label == 'be_z':
        return 'orange'
    elif label == 'be_x':
        return 'green'
    elif label == 'ke_y':
        return 'red'
    elif label == 'ke_z':
        return 'purple'
    elif label == 'ke_x':
        return 'brown'
    elif label == r"$\langle\frac{1}{2\tau}|\mathbf{u}|^2\rangle$":
        return 'black'
    elif label == 'jy':
        return 'green'
    elif label == 'jz':
        return 'purple'
    else:
        return np.random.rand(3,) # not sure if this works lol

def get_title():
    return r"Ro$=$" + str(Ro) + r"; $Pm=$" + str(Pm) + r"; $\nu=$" + str(nu)

def get_times(dir_list):
    times_lst = []
    for dir in dir_list:
        times_lst.append(0)
        files = glob.glob("{}scalars/*.h5".format(dir))
        last_index = len(files)
        for file in files:
            # if "_s{}.h5".format(last_index) in file:
            #     print("success")
            #     continue
            with h5py.File(file, "r") as f:
                sim_times = f['scales']['sim_time'][()]
                # print(sim_times)
                if max(sim_times) > times_lst[-1]:
                    times_lst[-1] = max(sim_times)
    return times_lst

def plot_energies(dir_list, time_list, prefix=''):
    initial_time = 0
    last_dir_indy = len(dir_list) - 1
    for indyboi in range(len(dir_list)):
        end_time = time_list[indyboi]
        dir = dir_list[indyboi]
        labelled = indyboi == last_dir_indy

        files = glob.glob("{}scalars/*.h5".format(dir))
        last_index = len(files)
        
        for file in files:
            # if "_s{}.h5".format(last_index) in file:
            #     print("success")
            #     continue

            with h5py.File(file, "r") as f:

                sim_times = f['scales']['sim_time'][()]

                be_y = f['tasks']['be_y'][()].ravel()
                be_z = f['tasks']['be_z'][()].ravel()
                be_x = f['tasks']['be_x'][()].ravel()
                ke_y = f['tasks']['ke_y'][()].ravel()
                ke_z = f['tasks']['ke_z'][()].ravel()
                ke_x = f['tasks']['ke_x'][()].ravel()

                # damp_power = -f['tasks']['damp_power'][()].ravel()

            labels = {
                'be_y' : be_y,
                'be_z' : be_z,
                'be_x' : be_x,
                'ke_y' : ke_y,
                'ke_z' : ke_z,
                'ke_x' : ke_x
                # r"$\langle\frac{1}{2\tau}|\mathbf{u}|^2\rangle$" : damp_power
            }

            
            # print('working')
            doPlot = True
            if max(sim_times) > end_time:
                # print(len(sim_times))
                # if 0 == next(x[0] for x in enumerate(sim_times) if x[1] > end_time):
                doPlot = False
            # sys.exit()


            for label in labels:
                if doPlot and prefix in label:
                    if labelled:
                        if prefix == '' and 'ke' in label:
                            plt.plot(initial_time + sim_times, labels[label], label=label, linestyle='dashed', color=get_color(label))
                        else:
                            plt.plot(initial_time + sim_times, labels[label], label=label, color=get_color(label))
                    else:
                        if prefix == '' and 'ke' in label:
                            plt.plot(initial_time + sim_times, labels[label], linestyle='dashed', color=get_color(label))
                        else:
                            plt.plot(initial_time + sim_times, labels[label], color=get_color(label))

        initial_time = end_time

    plt.legend()
    plt.title(get_title())
    # plt.title(r"$\tau=$" + str(tau))
    plt.xlabel("time")
    # plt.ylim(1e-3, 2.0)
    plt.ylabel("energy")
    plt.yscale('log')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.savefig('{}energies_{}.png'.format(dir, prefix))
    print('{}energies_{}.png'.format(dir, prefix))
    plt.close()

def plot_j(dir, offset=0):
    files = glob.glob("{}scalars/*.h5".format(dir))
    last_index = len(files)
    for file in files:
        if "_s{}.h5".format(last_index) in file:
            print("success")
            continue
        with h5py.File(file, "r") as f:

            sim_times = f['scales']['sim_time'][()]

            j = f['tasks']['j'][()]
            jy = j[:, 0, 0, 0, 0]
            jz = j[:, 1, 0, 0, 0]
            jx = j[:, 2, 0, 0, 0]

        labels = {
            # 'jx' : jx,
            'jy' : jy,
            'jz' : jz
        }


        for label in labels:

            if file == files[0]:
                plt.plot(offset + sim_times, labels[label], label=label, color=get_color(label))
            else:
                plt.plot(offset + sim_times, labels[label], color=get_color(label))


    plt.legend()
    plt.title(get_title())
    plt.xlabel("time")
    plt.ylabel("current density")
    # plt.yscale('log')
    figname = "{}{}j_mean.png".format(dir)
    plt.savefig(figname)
    print(figname)
    plt.close()

def plot_bdiff(dir, offset=0):
    files = glob.glob("{}scalars/*.h5".format(dir))
    last_index = len(files)
    for file in files:
        if "_s{}.h5".format(last_index) in file:
            print("success")
            continue
        with h5py.File(file, "r") as f:
            sim_times = f['scales']['sim_time'][()]
            bdiff = f['tasks']['be_y-be_z'][()].ravel()

        plt.plot(offset + sim_times, bdiff, color='purple')

    # plt.legend()
    plt.title(get_title())
    plt.xlabel("time")
    plt.ylabel("be_y - be_z")
    # plt.yscale('log')
    figname = "{}bdiff.png".format(dir)
    plt.savefig(figname)
    print(figname)
    plt.close()

from configparser import ConfigParser
try:
    filename = glob.glob("{}/*.cfg".format(sys.argv[1]))[0]
    config = ConfigParser()
    config.read(str(filename))
    Pm = config.getfloat("parameters", "Pm")
    Ro = config.getfloat("parameters", "Ro")
    tau = abs(config.getfloat("parameters", "tau"))
    nu = config.getfloat("parameters", "nu")
except:
    print("failed to read config file. Please supply run suffix")
    raise

# parent1 = sys.argv[1] + '/data/'
# parent2 = sys.argv[2] + '/data/'

parents = ["{}/data/1/".format(ss) for ss in sys.argv[1:]]
raw_times = get_times(parents)
rounded_times = [50*int(raw_time/50) for raw_time in raw_times]

plot_energies(parents, rounded_times, prefix=sys.argv[-1])

size = CW.size
rank = CW.rank 
targets = []
seeds = []

N = len(list(range(rank, len(seeds), size)))
