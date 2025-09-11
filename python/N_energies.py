import numpy as np
# import dedalus.public as d3
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from mpi4py import MPI
CW = MPI.COMM_WORLD

def get_label(arg):
    if arg == 'jy':
        return r"$\langle j_y \rangle$"
    elif arg == 'jz':
        return r"$\langle j_z \rangle$"
    elif arg == 'be_y':
        return r"$0.5\langle |b_y^2| \rangle$"
    elif arg == 'be_z':
        return r"$0.5\langle |b_z^2| \rangle$"
    elif arg == 'be_x':
        return r"$0.5\langle |b_x^2| \rangle$"
    elif arg == 'ke_y':
        return r"$0.5\langle |u_y^2| \rangle$"
    elif arg == 'ke_z':
        return r"$0.5\langle |u_z^2| \rangle$"
    elif arg == 'ke_x':
        return r"$0.5\langle |u_x^2| \rangle$"


    else:
        raise

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
    return None
    # return "title under construction..."
    # return r"Ro$=$" + str(Ro) + r"; $Pm=$" + str(Pm) + r"; $\nu=$" + str(nu)

suffixes = sys.argv[1:]
print("suffixes = {}".format(suffixes))
# sys.exit()
from configparser import ConfigParser
try:
    filename = glob.glob("{}/*.cfg".format(sys.argv[1]))[0]
    config = ConfigParser()
    config.read(str(filename))
    # Pm = config.getfloat("parameters", "Pm")
    # nu = config.getfloat("parameters", "nu")
    # try:
    #     Ro = config.getfloat("parameters", "Ro")
    # except:
    #     t_vec = eval(eval(config.get("parameters", "t_vec")))
    #     Ro_vec = (eval(config.get("parameters", "Ro_vec")))
    #     Ro = (Ro_vec[0], Ro_vec[-1])
    # icContinuation = config.getfloat("parameters", "Ro")

except:
    print("failed to read config file. Please supply run suffix")
    raise

prefix = ''
parents = [suf + '/' for suf in suffixes]
# plot_energies(parents, prefix='')

init_t = 0.0
end_sim_times = []
for dir_name in parents:
    files = glob.glob("{}checkpoint/*.h5".format(dir_name))
    last_index = len(files)
    for file in files:
        if 'checkpoint_s{}.h5'.format(last_index) in file:
            with h5py.File(file, "r") as f:
                end_sim_times.append(f['scales']['sim_time'][()][0])
try:
    end_sim_times[-1] = 0.0
except:
    end_sim_times.append(0.0)
# print(end_sim_times)
# sys.exit()
delta_t = init_t = 0.0
for indy, dir_name in enumerate(parents):
    # print(end_sim_times[indy])
    print(init_t)
    files = glob.glob("{}scalars/*.h5".format(dir_name))
    last_index = len(files)
    for file in files:
            # print("success")
            # continue

        with h5py.File(file, "r") as f:

            sim_times = init_t + f['scales']['sim_time'][()]
            if sim_times[0] >= end_sim_times[indy] and indy != len(parents) - 1:
                continue

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

        start_ind = 0
        passed_checkpoint = not np.any(sim_times < end_sim_times[indy])
        for label in labels:

            if file == files[0] and init_t == 0:
                if prefix == '' and 'ke' in label:
                    plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], label=get_label(label), linestyle='dashed', color=get_color(label))
                else:
                    plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], label=get_label(label), color=get_color(label))
            # elif passed_checkpoint:
            #     if prefix == '' and 'ke' in label:
            #         plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], linestyle='dashed', color=get_color(label))
            #     else:
            #         plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], color=get_color(label))

            else:
                if prefix == '' and 'ke' in label:
                    plt.plot(sim_times[start_ind:], labels[label][start_ind:], linestyle='dashed', color=get_color(label))
                else:
                    plt.plot(sim_times[start_ind:], labels[label][start_ind:], color=get_color(label))

        # if 'scalars_s{}.h5'.format(last_index) in file:
        if abs(sim_times[-1] - end_sim_times[indy]) <= 0.05:
            try:
                # print(sim_times)
                delta_t = sim_times[-1]
            except:
                delta_t = 0.0
    init_t += delta_t
    # if indy == 0:
    #     break

plt.legend(framealpha=1.0, ncol=2, loc='center left', bbox_to_anchor=(1, 0.5))
plt.title(get_title())
# plt.title(r"$\tau=$" + str(tau))
plt.xlabel("time")
# plt.ylim(1e-3, 2.0)
plt.ylabel("energy")
plt.yscale('log')
for dir_name in parents:
    plt.savefig('{}combined_{}.png'.format(dir_name, prefix))
    print('{}combined_{}.png'.format(dir_name, prefix))
plt.close()
