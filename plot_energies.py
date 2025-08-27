import numpy as np
# import dedalus.public as d3
import matplotlib.pyplot as plt
import h5py
import glob
import sys

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

def plot_energies(dir, prefix='ke'):
    files = glob.glob("{}scalars/*.h5".format(dir))
    last_index = len(files)
    for file in files:
        if "_s{}.h5".format(last_index) in file:
            print("success")
            # continue

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

        start_ind = 0
        for label in labels:
            if prefix in label:
                if file == files[0]:
                    if prefix == '' and 'ke' in label:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], label=get_label(label), linestyle='dashed', color=get_color(label))
                    else:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], label=get_label(label), color=get_color(label))
                else:
                    if prefix == '' and 'ke' in label:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], linestyle='dashed', color=get_color(label))
                    else:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], color=get_color(label))
    plt.legend(framealpha=1.0, ncol=3)
    plt.title(get_title())
    # plt.title(r"$\tau=$" + str(tau))
    plt.xlabel("time")
    # plt.ylim(1e-3, 2.0)
    plt.ylabel("energy")
    plt.yscale('log')
    plt.savefig('{}energies_{}.png'.format(dir, prefix))
    print('{}energies_{}.png'.format(dir, prefix))
    plt.close()

def plot_2energies(dir, prefix='ke'):
    plt.figure(figsize=(4, 3))
    files = glob.glob("{}scalars/*.h5".format(dir))
    last_index = len(files)
    for file in files:
        if "_s{}.h5".format(last_index) in file:
            print("success")
            # continue

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
            'ke_y' : ke_y,
            'be_y' : be_y            # r"$\langle\frac{1}{2\tau}|\mathbf{u}|^2\rangle$" : damp_power
        }

        start_ind = 0
        for label in labels:
            if prefix in label:
                if file == files[0]:
                    if prefix == '' and 'ke' in label:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], label=get_label(label), linestyle='solid', color=get_color(label))
                    else:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], label=get_label(label), color=get_color(label))
                else:
                    if prefix == '' and 'ke' in label:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], linestyle='solid', color=get_color(label))
                    else:
                        plt.plot(sim_times[start_ind:], labels[label][start_ind:], color=get_color(label))
    plt.legend(framealpha=0.0)
    # plt.legend(framealpha=0.0, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 0.9))
    plt.title(get_title())

    # plt.title(r"$\tau=$" + str(tau))
    plt.xlabel("time")
    # plt.xlim(1000, 1500)
    # plt.ylim(3e-2, 0.15)
    # plt.ylabel("energy")
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig('{}energies2_{}.png'.format(dir, prefix))
    print('{}energies2_{}.png'.format(dir, prefix))
    plt.close()

def plot_j(dir):
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


        for qb, label in enumerate(labels):
            linethickness = 0.3
            if qb == 0:
                if file == files[0]:
                    plt.plot(sim_times, labels[label], label=get_label(label), color=get_color(label), linewidth=linethickness)
                else:
                    plt.plot(sim_times, labels[label], color=get_color(label), linewidth=linethickness)
            else:
                linestyle = 'dashed'
                if file == files[0]:
                    plt.plot(sim_times, labels[label], label=get_label(label), color=get_color(label), linestyle=linestyle, linewidth=linethickness)
                else:
                    plt.plot(sim_times, labels[label], color=get_color(label), linestyle=linestyle, linewidth=linethickness)

    plt.legend()
    plt.title(get_title())
    plt.xlabel("time")
    plt.ylabel("current density")
    # plt.yscale('log')
    figname = "{}j_mean.png".format(dir)
    plt.savefig(figname)
    print(figname)
    plt.close()

def plot_task(dir, task, logscale=False):
    files = glob.glob("{}scalars/*.h5".format(dir))
    files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    last_index = len(files)
    times_data = []
    data_data = []
    for file in files:
        # print(file)
        if "_s{}.h5".format(last_index) in file:
            print("success")
            continue
        with h5py.File(file, "r") as f:
            times_data += f['scales']['sim_time'][()].tolist()
            data_data += f['tasks'][task][()].ravel().tolist()

    if 'proj' in task:
        plt.plot(times_data, 1 - np.array(data_data), color='purple')
        plt.ylim(5e-3, 3)
    else:
        plt.plot(times_data, data_data, color='purple')

    # plt.legend()
    plt.title(get_title())
    plt.xlabel("time")
    plt.ylabel(task)
    if logscale:
        plt.yscale('log')
    if task == "udiff":
        plt.ylim(1e-16, 1e-1)
    elif task == 'ke_y':
        plt.ylim(7e-2, 8e-2)
    # plt.yscale('log')
    plt.tight_layout()
    figname = "{}{}.png".format(dir,task)
    plt.savefig(figname)
    print(figname)
    plt.close()

def plot_tasks(dir, tasks, logscale=False):
    files = glob.glob("{}scalars/*.h5".format(dir))
    files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    last_index = len(files)
    times_data = []
    data_dict = dict()
    for task in tasks:
        data_dict[task] = []
    for file in files:
        # print(file)
        if "_s{}.h5".format(last_index) in file:
            print("success")
            continue
        with h5py.File(file, "r") as f:
            times_data += f['scales']['sim_time'][()].tolist()
            for task in tasks:
                data_dict[task] += f['tasks'][task][()].ravel().tolist()

    for task in tasks:
        plt.plot(times_data, data_dict[task], label=task)

    plt.legend()
    plt.title(get_title())
    plt.xlabel("time")
    plt.ylabel(task)
    if logscale:
        plt.yscale('log')
    if task == "udiff":
        plt.ylim(1e-16, 1e-1)

    # plt.yscale('log')
    figname = "{}{}.png".format(dir,task)
    plt.savefig(figname)
    print(figname)
    plt.close()


def plot_23_phase(dir, prefix='ke'):
    plt.figure(figsize=(4, 3))

    files = glob.glob("{}scalars/*.h5".format(dir))
    last_index = len(files)
    for file in files:
        if "_s{}.h5".format(last_index) in file:
            print("success")
            # continue
        labels = {}
        with h5py.File(file, "r") as f:
            Nmodes = 4
            ke_modes = []
            for ki in range(Nmodes):
                ke_modes.append(f['tasks']['{}_mode{}'.format(prefix, ki + 1)][()].ravel())
                labels['{}_mode{}'.format(prefix, ki + 1)] = ke_modes[-1]

        start_ind = 0
        plt.plot(labels['ke_mode2'], labels['ke_mode3'], color='magenta')
    # plt.legend(framealpha=0.0, loc='best')
    plt.title(get_title())
    # plt.title(r"$\tau=$" + str(tau))
    plt.xlabel("mode 2")
    # plt.ylim(5e-3, 5e-2)
    plt.ylabel("mode 3")
    # plt.yscale('log')
    plt.tight_layout()
    plt.savefig('{}{}phase_23.png'.format(dir, prefix))
    print('{}{}phase_23.png'.format(dir, prefix))
    plt.close()

def plot_zmodes(dir, prefix='ke'):
    plt.figure(figsize=(4, 3))
    def color10(label):
        ki = int(label[7:])
        if ki == 1:
            return 'black'
        elif ki == 2:
            return 'grey'
        elif ki == 3:
            return 'lime'
        elif ki == 4:
            return  'pink'
        elif ki == 5:
            return 'purple'
        elif ki == 6:
            return 'navy'
        elif ki == 7:
            return 'green'
        elif ki == 8:
            return 'orange'
        elif ki == 9:
            return 'yellow'
        elif ki == 10:
            return 'brown'
        else:
            raise
    files = glob.glob("{}scalars/*.h5".format(dir))
    last_index = len(files)
    for file in files:
        if "_s{}.h5".format(last_index) in file:
            print("success")
            # continue

        labels = {}
        with h5py.File(file, "r") as f:

            sim_times = f['scales']['sim_time'][()]

            Nmodes = 4
            ke_modes = []
            for ki in range(Nmodes):
                ke_modes.append(f['tasks']['{}_mode{}'.format(prefix, ki + 1)][()].ravel())
                labels['{}_mode{}'.format(prefix, ki + 1)] = ke_modes[-1]

            # damp_power = -f['tasks']['damp_power'][()].ravel()

            # 'be_y' : be_y,
            # 'be_z' : be_z,
            # 'be_x' : be_x,
            # 'ke_y' : ke_y,
            # 'ke_z' : ke_z,
            # 'ke_x' : ke_x
            # r"$\langle\frac{1}{2\tau}|\mathbf{u}|^2\rangle$" : damp_power
        # }

        start_ind = 0
        for label in labels:
            if file == files[0]:
                if 'ke' in label:
                    plt.plot(sim_times[start_ind:], labels[label][start_ind:], label=label, linestyle='solid', color=color10(label))
                else:
                    plt.plot(sim_times[start_ind:], labels[label][start_ind:], label=label, color=color10(label))
            else:
                if 'ke' in label:
                    plt.plot(sim_times[start_ind:], labels[label][start_ind:], linestyle='solid', color=color10(label))
                else:
                    plt.plot(sim_times[start_ind:], labels[label][start_ind:], color=color10(label))
    plt.legend(framealpha=0.0, loc='best')
    plt.gca().set_ylim(bottom=0)
    plt.title(get_title())
    # plt.title(r"$\tau=$" + str(tau))
    plt.xlabel("time")
    # plt.ylim(5e-3, 5e-2)
    # plt.ylabel("energy")
    # plt.yscale('log')
    plt.tight_layout()
    plt.savefig('{}{}zmodes.png'.format(dir, prefix))
    print('{}{}zmodes.png'.format(dir, prefix))
    plt.close()

def plot_invariant(dir):
    files = glob.glob("{}scalars/*.h5".format(dir))
    last_index = len(files)
    for file in files:
        if "_s{}.h5".format(last_index) in file:
            print("success")
            continue

        labels = {}
        with h5py.File(file, "r") as f:

            sim_times = f['scales']['sim_time'][()]
            ke_mode2 = f['tasks']['ke_mode2'][()].ravel()
            ke_mode3 = f['tasks']['ke_mode3'][()].ravel()
            ke_y = f['tasks']['ke_y'][()].ravel()
            ke_z = f['tasks']['ke_z'][()].ravel()
            ke_x = f['tasks']['ke_x'][()].ravel()
            be_y = f['tasks']['be_y'][()].ravel()
            be_z = f['tasks']['be_z'][()].ravel()
            be_x = f['tasks']['be_x'][()].ravel()
            invariant = (be_y + be_z + be_x) + (ke_y + ke_z + ke_x)
            # invariant = be_y * ke_x**0.5
            # invariant = (be_y) / (ke_mode3 + 0.1*ke_mode2)

        start_ind = 0
        plt.plot(sim_times[start_ind:], invariant[start_ind:], color='black')

    plt.title(get_title())
    # plt.title(r"$\tau=$" + str(tau))
    plt.xlabel("time")
    # plt.ylim(5e-3, 5e-2)
    plt.ylabel("be_y / ke_3")
    # plt.yscale('log')
    plt.savefig('{}invar.png'.format(dir))
    print('{}invar.png'.format(dir))
    plt.close()



from configparser import ConfigParser
try:
    filename = glob.glob("{}/*.cfg".format(sys.argv[1]))[0]
    config = ConfigParser()
    config.read(str(filename))
    # Pm = config.getfloat("parameters", "Pm")
    # try:
    #     Ro = config.getfloat("parameters", "Ro")
    # except:
    #     t_vec = eval(eval(config.get("parameters", "t_vec")))
    #     Ro_vec = (eval(config.get("parameters", "Ro_vec")))
    #     Ro = (Ro_vec[0], Ro_vec[-1])
    # icContinuation = config.getfloat("parameters", "Ro")

    # nu = config.getfloat("parameters", "nu")
except:
    print("failed to read config file. Please supply run suffix")
    raise

parent = sys.argv[1]
targets = []
seeds = []
# for folder in glob.glob(parent + "*/"):
#     folder = folder.replace(parent, "")
#     folder = folder.replace("/", "")
#     if folder.isdigit():
#         targets.append(folder)
#         seeds.append(int(folder))

# N = len(list(range(rank, len(seeds), size)))

# for i in range(rank, len(seeds), size):
# if (rank == 0):
    # print('plotting from simulation {} / {} on rank 0...'.format(i, N))
# try:
#     plot_tasks(parent + '/', ['phase_ux', 'phase_decoy'], logscale=False)
# except Exception as e:
#     print(e)
#     print('failed phases')
    
# try:
#     plot_task(parent + '/', "B0_energy_ratio_yz", logscale=False)
# except Exception as e:
#     print(e)
#     print('failed B0_energy_ratio_yz')

# try:
#     plot_task(parent + '/', "B0_energy_ratio_y", logscale=False)
# except Exception as e:
#     print(e)
#     print('failed B0_energy_ratio_y')

# try:
#     plot_task(parent + '/', 'projA', logscale=True)
# except Exception as e:
#     print('failed to plot projection projA')

# try:
#     plot_task(parent + '/', 'projb', logscale=True)
# except Exception as e:
#     print('failed to plot projection projb')

try:
    plot_energies(parent + '/', prefix='')
except Exception as e:
    print('failed to plot energies')
try:
    plot_2energies(parent + '/', prefix='')
except Exception as e:
    print('failed to plot energies')
    print(e)
try:
    plot_invariant(parent + '/')
except Exception as e:
    print('failed to plot invariant')
try:
    plot_zmodes(parent + '/')
except Exception as e:
    print(e)
    print('failed to plot zmodes')

try:
    plot_23_phase(parent + '/')
except Exception as e:
    print(e)
    print('failed to plot phase_23')

try:
    plot_zmodes(parent + '/', prefix='be')
except Exception as e:
    print(e)
    print('failed to plot zmodes')

try:
    plot_task(parent + '/', 'keff', logscale=False)
except Exception as e:
    print('failed to plot keff')

try:
    plot_task(parent + '/', 'ke_y', logscale=False)
except Exception as e:
    print('failed to plot ke_y')

try:
    plot_task(parent + '/', 'bhelicity', logscale=False)
except Exception as e:
    print('failed to plot bhelicity')

try:
    plot_task(parent + '/', 'be', logscale=True)
except Exception as e:
    print('failed to plot be')

try:
    plot_task(parent + '/', 'be_y', logscale=True)
except Exception as e:
    print('failed to plot be_y')

try:
    plot_task(parent + '/', 'ke_x', logscale=True)
except Exception as e:
    print('failed to plot ke_x')

try:
    plot_task(parent + '/', 'udiff', logscale=True)
except Exception as e:
    print('failed to plot udiff')

try:
    plot_task(parent + '/', 'Adiff', logscale=True)
except Exception as e:
    print('failed to plot udiff')

# sys.exit()
# plot_j(parent + '/')
# plot_bdiff(parent + '/')
try:
    plot_task(parent + '/', 'omega', logscale=False)
except Exception as e:
    print('failed to plot omega')

try:
    plot_task(parent + '/', 'omega', logscale=False)
except Exception as e:
    print('failed to plot omega')

try:
    plot_task(parent + '/', 'proj_A0')
except Exception as e:
    print('failed to plot projection A0')

try:
    plot_task(parent + '/', 'proj_norm_A0')
except Exception as e:
    print('failed to plot projection A0 normalized')

try:
    plot_task(parent + '/', 'be_y-be_z')
except Exception as e:
    print('failed to plot b-diff')
    
try:
    plot_task(parent + '/', 'Rossby')
except Exception as e:
    print('failed to plot Rossby')
    
    
try:
    plot_task(parent + '/', 'Rayleigh')
except Exception as e:
    print('failed to plot Rayleigh')
    
    
    
try:
    plot_task(parent + '/', 'AdotA_mean')
except Exception as e:
    print('failed to plot AdotA_mean')
    