import numpy as np
# import dedalus.public as d3
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from mpi4py import MPI
CW = MPI.COMM_WORLD
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 6
mpl.rcParams['figure.figsize'] = (6, 4)


justOne = True
doTruncate = True
truncateSecond = False
truncate_time = 60
# suffixes = suffixes1 = ['nl_Ro5_nu0p0076923', 'k1_twoRo5_nu0p0076923_cp7']
suffixes = suffixes1 = ['k3_two_stable_Ro5_nu0p0076923']
# suffixes = suffixes1 = ['k3_three_Ro5_nu0p0076923']
# suffixes = suffixes1 = ['B0p005', 'k1_B0p005']
# suffixes2 = ['case1_noise', 'k1_case1three']
delta_t = 0
# truncate_time = 80
# suffixes = suffixes1 = ['nl_Ro5_nu0p0076923', 'k1_twoRo5_nu0p0076923_cp9']

def get_color(label):
    if label == 0:
        return 'cyan'
    elif label == 1:
        return 'magenta'

def get_title():
    return None
    # return "title under construction..."
    # return r"Ro$=$" + str(Ro) + r"; $Pm=$" + str(Pm) + r"; $\nu=$" + str(nu)

# suffixes = sys.argv[1:]
print("suffixes = {}".format(suffixes))

# sys.exit()
from configparser import ConfigParser
try:
    filename = glob.glob("{}/*.cfg".format(suffixes[0]))[0]
    config = ConfigParser()
    config.read(str(filename))
    Pm = config.getfloat("parameters", "Pm")
    nu = config.getfloat("parameters", "nu")

except:
    print("failed to read config file. Please supply run suffix")
    raise








prefix = ''
parents = [suf + '/data/1/' for suf in suffixes]
# plot_energies(parents, prefix='')

init_init_t = init_t = 0.0
end_sim_times = [truncate_time]
for dir_name in parents:
    files = glob.glob("{}checkpoint/*.h5".format(dir_name))
    last_index = len(files)
    for file in files:
        if 'checkpoint_s{}.h5'.format(last_index) in file:
            with h5py.File(file, "r") as f:
                end_sim_times.append(f['scales']['sim_time'][()][0])
end_sim_times[-1] = 0.0

# print(end_sim_times)
# sys.exit()
init_t = 0.0
for indy, dir_name in enumerate(parents):
    # print(end_sim_times[indy])
    print(init_t)
    init_init_t = init_t
    files = glob.glob("{}scalars/*.h5".format(dir_name))
    last_index = len(files)
    for file in files:
            # print("success")
            # continue

        with h5py.File(file, "r") as f:

            sim_times = init_t + f['scales']['sim_time'][()]
            if sim_times[0] >= end_sim_times[indy] and indy != len(end_sim_times) - 1:
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
            # 'be_z' : be_z,
            # 'be_x' : be_x,
            # 'ke_y' : ke_y,
            # 'ke_z' : ke_z,
            # 'ke_x' : ke_x
            # r"$\langle\frac{1}{2\tau}|\mathbf{u}|^2\rangle$" : damp_power
        }

        start_ind = 0
        passed_checkpoint = not np.any(sim_times < end_sim_times[indy])
        for label in labels:

            if file == files[0] and init_t == 0:
                if prefix == '' and 'ke' in label:
                    plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], linestyle='dashed', color=get_color(0))
                else:
                    plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], color=get_color(0))
            # elif passed_checkpoint:
            #     if prefix == '' and 'ke' in label:
            #         plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], linestyle='dashed', color=get_color(label))
            #     else:
            #         plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], color=get_color(label))

            else:
                if prefix == '' and 'ke' in label:
                    plt.plot(sim_times[start_ind:], labels[label][start_ind:], linestyle='dashed', color=get_color(0))
                else:
                    plt.plot(sim_times[start_ind:], labels[label][start_ind:], color=get_color(0))

        # if 'scalars_s{}.h5'.format(last_index) in file:
        if abs(sim_times[-1] - end_sim_times[indy]) <= 0.05:
            try:
                # print(sim_times)
                delta_t = sim_times[-1]
            except:
                delta_t = 0.0
    init_t += delta_t
    if init_init_t > 0 and not doTruncate:
        print('yes')
        plt.axvline(x=init_init_t, color='grey', linestyle='dashed')
    # if indy == 0:
    #     break

if not justOne:
    suffixes = suffixes2



    prefix = ''
    parents = [suf + '/data/1/' for suf in suffixes]
    # plot_energies(parents, prefix='')

    init_init_t = init_t = 0.0
    end_sim_times = []
    for dir_name in parents:
        files = glob.glob("{}checkpoint/*.h5".format(dir_name))
        last_index = len(files)
        for file in files:
            if 'checkpoint_s{}.h5'.format(last_index) in file:
                with h5py.File(file, "r") as f:
                    end_sim_times.append(f['scales']['sim_time'][()][0])
    end_sim_times[-1] = 0.0

    # print(end_sim_times)
    # sys.exit()
    init_t = 0.0
    for indy, dir_name in enumerate(parents):
        # print(end_sim_times[indy])
        print(init_t)
        init_init_t = init_t
        files = glob.glob("{}scalars/*.h5".format(dir_name))
        last_index = len(files)
        for file in files:
                # print("success")
                # continue

            with h5py.File(file, "r") as f:

                sim_times = init_t + f['scales']['sim_time'][()]
                if sim_times[0] >= end_sim_times[indy] and indy != len(parents) - 1:
                    continue
                if truncateSecond and sim_times[0] >= truncate_time:
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
                # 'be_z' : be_z,
                # 'be_x' : be_x,
                # 'ke_y' : ke_y,
                # 'ke_z' : ke_z,
                # 'ke_x' : ke_x
                # r"$\langle\frac{1}{2\tau}|\mathbf{u}|^2\rangle$" : damp_power
            }

            start_ind = 0
            passed_checkpoint = not np.any(sim_times < end_sim_times[indy])
            for label in labels:

                if file == files[0] and init_t == 0:
                    if prefix == '' and 'ke' in label:
                        plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], linestyle='dashed', color=get_color(1))
                    else:
                        plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], color=get_color(1))
                # elif passed_checkpoint:
                #     if prefix == '' and 'ke' in label:
                #         plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], linestyle='dashed', color=get_color(label))
                #     else:
                #         plt.plot(sim_times[sim_times > end_sim_times[indy]], labels[label][sim_times > end_sim_times[indy]], color=get_color(label))

                else:
                    if prefix == '' and 'ke' in label:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], linestyle='dashed', color=get_color(1))
                    else:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], color=get_color(1))

            # if 'scalars_s{}.h5'.format(last_index) in file:
            if abs(sim_times[-1] - end_sim_times[indy]) <= 0.05:
                try:
                    # print(sim_times)
                    delta_t = sim_times[-1]
                except:
                    delta_t = 0.0
        init_t += delta_t
        if init_init_t > 0 and not doTruncate:
            print('yes')
            plt.axvline(x=init_init_t, color='grey', linestyle='dashed')
        # if indy == 0:
        #     break

# plt.legend(framealpha=1.0, ncol=3)
plt.title("azimuthal magnetic energy")
# plt.title(r"$\tau=$" + str(tau))
plt.xlabel("time")
# plt.xlim(0, 100.0)
if doTruncate:
    plt.xlim(0, truncate_time)
# plt.ylim(1e-3, 2.0)
# plt.ylabel("energy")
plt.yscale('log')
for dir_name in parents:
    if justOne:
        figname = '{}Oneclen_{}.png'.format(dir_name, prefix)
    else:
        figname = '{}Bothclen_{}.png'.format(dir_name, prefix)
    if doTruncate:
        figname = figname[:-4] + 'trunc.png'
    plt.savefig(figname, dpi=300)
    print(figname)
plt.close()
