import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
import os
import sys
from configparser import ConfigParser

# Combines the scalar-diagnostic plots from plot_energies.py across a chain of
# restarted runs (e.g. `Ro3p5_Lz8pi_Ny16_Rm175` -> `Ro3p5_Lz8pi_Ny16_Rm175_RSTRT1`).
#
# Every restart in this codebase reloads the field state from a checkpoint but
# resets solver.sim_time to zero (see rpcf-mhd.py), so each run directory's own
# scalar files are timestamped from 0 rather than continuing the true sim time.
# To combine them into one continuous timeline:
#   - every run except the last is trimmed to the sim_time of its own final
#     checkpoint (the point the next run in the chain actually restarted from,
#     since anything logged after that checkpoint was never resumed from)
#   - each subsequent run's local sim_time is shifted by the cumulative trim
#     time of the runs before it
#
# Usage:
#   python3 plot_energies_combined.py <dir1> <dir2> [<dir3> ...]
# e.g.
#   python3 plot_energies_combined.py Ro3p5_Lz8pi_Ny16_Rm175 Ro3p5_Lz8pi_Ny16_Rm175_RSTRT1
#
# Figures are saved under the same names as plot_energies.py with "_combined"
# inserted before the extension, written into every directory passed in.

def get_label(arg):
    if arg == 'jy':
        return r"$\langle j_y \rangle$"
    elif arg == 'jz':
        return r"$\langle j_z \rangle$"
    elif arg == 'be_y':
        return r"$0.5\langle b_y^2 \rangle$"
    elif arg == 'be_z':
        return r"$0.5\langle b_z^2 \rangle$"
    elif arg == 'be_x':
        return r"$0.5\langle b_x^2 \rangle$"
    elif arg == 'ke_y':
        return r"$0.5\langle u_y^2 \rangle$"
    elif arg == 'ke_z':
        return r"$0.5\langle u_z^2 \rangle$"
    elif arg == 'ke_x':
        return r"$0.5\langle u_x^2 \rangle$"
    else:
        raise

def get_color(label):
    if label == 'be_y':
        return 'deepskyblue'
    elif label == 'be_z':
        return 'orange'
    elif label == 'be_x':
        return 'green'
    elif label == 'ke_y':
        return 'magenta'
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
        return np.random.rand(3,)

def get_title():
    return None


# ---------------------------------------------------------------------------
# multi-run data loading
# ---------------------------------------------------------------------------

def get_checkpoint_final_time(dir):
    """sim_time of the last (highest-index) checkpoint saved in `dir`."""
    files = glob.glob("{}checkpoint/checkpoint_s*.h5".format(dir))
    if not files:
        return None
    def idx(f):
        return int(os.path.basename(f).split('checkpoint_s')[-1].split('.h5')[0])
    files.sort(key=idx)
    with h5py.File(files[-1], "r") as f:
        return float(f['scales']['sim_time'][()][0])

def sorted_scalar_files(dir):
    files = glob.glob("{}scalars/*.h5".format(dir))
    files.sort(key=lambda f: int(''.join(filter(str.isdigit, os.path.basename(f)))))
    return files

def load_series(dirs, tasks):
    """
    Load and concatenate the requested scalar tasks (plus sim_time) across a
    chain of run directories, trimming every run but the last to its own
    final checkpoint time and offsetting later runs' local sim_time so the
    result reads as one continuous timeline.
    """
    series = {t: [] for t in tasks}
    series['sim_time'] = []
    offset = 0.0
    for i, dir in enumerate(dirs):
        trim_time = None
        if i < len(dirs) - 1:
            trim_time = get_checkpoint_final_time(dir)
            if trim_time is None:
                print("warning: no checkpoints found in {}, not trimming".format(dir))
        for file in sorted_scalar_files(dir):
            with h5py.File(file, "r") as f:
                sim_times = f['scales']['sim_time'][()]
                mask = np.ones(sim_times.shape, dtype=bool)
                if trim_time is not None:
                    mask &= sim_times <= trim_time
                if not mask.any():
                    continue
                series['sim_time'].append(sim_times[mask] + offset)
                for t in tasks:
                    series[t].append(f['tasks'][t][()].ravel()[mask])
        if trim_time is not None:
            offset += trim_time
    for k in series:
        series[k] = np.concatenate(series[k]) if series[k] else np.array([])
    return series

def combined_name(filename, tag='_combined'):
    base, ext = filename.rsplit('.', 1)
    return "{}{}.{}".format(base, tag, ext)

def savefig_combined(dirs, filename):
    fname = combined_name(filename)
    path = dirs[0] + fname
    plt.savefig(path)
    print(path)


# ---------------------------------------------------------------------------
# plots (mirrors plot_energies.py, operating on a list of dirs)
# ---------------------------------------------------------------------------

def plot_energies(dirs, prefix=''):
    plt.figure(figsize=(4, 3))
    tasks = ['be_y', 'be_z', 'be_x', 'ke_y', 'ke_z', 'ke_x']
    data = load_series(dirs, tasks)
    sim_times = data['sim_time']
    for label in tasks:
        if prefix in label:
            if prefix == '' and 'ke' in label:
                plt.plot(sim_times, data[label], label=get_label(label), linestyle='dashed', color=get_color(label))
            else:
                plt.plot(sim_times, data[label], label=get_label(label), color=get_color(label))
    plt.legend(framealpha=1.0, ncol=3)
    plt.title(get_title())
    plt.xlabel("time")
    plt.ylabel("energy")
    plt.yscale('log')
    savefig_combined(dirs, 'energies_{}.png'.format(prefix))
    plt.close()

def plot_2energies(dirs, prefix='ke', config=None):
    try:
        Rm = config.getfloat('parameters', 'Rm')
    except:
        Rm = None
    plt.figure(figsize=(4, 3))
    tasks = ['ke_y', 'be_y']
    data = load_series(dirs, tasks)
    sim_times = data['sim_time']
    max_sim_time = sim_times.max() if len(sim_times) else 0
    for label in tasks:
        if prefix in label:
            if prefix == '' and 'ke' in label:
                plt.plot(sim_times, data[label], label=get_label(label), linestyle='solid', color=get_color(label))
            else:
                plt.plot(sim_times, data[label], label=get_label(label), color=get_color(label))
    plt.legend(framealpha=0.0)
    if any('1e3' in d for d in dirs):
        plt.ylabel("nonlinear")
    plt.xlabel("time")
    plt.xlim(left=0)
    max_sim_time = min(max_sim_time, 1e4)
    plt.xlim(right=max_sim_time)
    if not any('seed' in d for d in dirs):
        plt.xlim(right=9999)
        plt.ylim(bottom=1e-6)
    else:
        plt.title(r"Rm = " + str(int(Rm)))
    plt.yscale('log')
    plt.tight_layout()
    savefig_combined(dirs, 'energies2_{}.png'.format(prefix))
    savefig_combined(dirs, 'energies2_{}.pdf'.format(prefix))
    plt.close()

def plot_invariant(dirs):
    tasks = ['ke_mode2', 'ke_mode3', 'ke_y', 'ke_z', 'ke_x', 'be_y', 'be_z', 'be_x']
    data = load_series(dirs, tasks)
    invariant = (data['be_y'] + data['be_z'] + data['be_x']) + (data['ke_y'] + data['ke_z'] + data['ke_x'])
    plt.plot(data['sim_time'], invariant, color='black')
    plt.title(get_title())
    plt.xlabel("time")
    plt.ylabel("be_y / ke_3")
    savefig_combined(dirs, 'invar.png')
    plt.close()

def plot_12_phase(dirs, prefix='ke'):
    plt.figure(figsize=(4, 3))
    tasks = ['{}_mode1'.format(prefix), '{}_mode2'.format(prefix)]
    data = load_series(dirs, tasks)
    plt.plot(data[tasks[0]], data[tasks[1]], color='magenta')
    plt.title(get_title())
    plt.xlabel("mode 1")
    plt.ylabel("mode 2")
    plt.tight_layout()
    savefig_combined(dirs, '{}phase_12.png'.format(prefix))
    plt.close()

def plot_zmodes(dirs, prefix='ke'):
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
            return 'pink'
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

    Nmodes = 10
    tasks = ['{}_mode{}'.format(prefix, ki + 1) for ki in range(Nmodes)]
    data = load_series(dirs, tasks)
    sim_times = data['sim_time']
    for label in tasks:
        if 'ke' in label:
            plt.plot(sim_times, data[label], label=label, linestyle='solid', color=color10(label))
        else:
            plt.plot(sim_times, data[label], label=label, color=color10(label))
    plt.legend(framealpha=0.0, loc='best')
    plt.gca().set_ylim(bottom=0)
    plt.title(get_title())
    plt.xlabel("time")
    plt.tight_layout()
    savefig_combined(dirs, '{}zmodes.png'.format(prefix))
    plt.close()

def plot_task(dirs, task, logscale=False):
    data = load_series(dirs, [task])
    times_data = data['sim_time']
    data_data = data[task]
    if 'proj' in task:
        plt.plot(times_data, 1 - np.array(data_data), color='purple')
        plt.ylim(5e-3, 3)
    else:
        plt.plot(times_data, data_data, color='purple')
    if 'Adiff' in task:
        plt.xlim(0, 23.33 * 40)
        plt.gca().set_ylim(bottom=1e-6)
    plt.title(get_title())
    plt.xlabel("time")
    plt.ylabel(task)
    if logscale:
        plt.yscale('log')
    if task == "udiff":
        plt.ylim(1e-16, 1e-1)
    elif task == 'ke_y':
        plt.ylim(7e-2, 8e-2)
    plt.tight_layout()
    savefig_combined(dirs, '{}.png'.format(task))
    plt.close()


def find_restart_chain(base_dir):
    """
    Auto-discover a restart chain from a base run directory by looking for
    <base_dir>_RSTRT1, <base_dir>_RSTRT2, ... (this project's restart naming
    convention), stopping at the first missing index.
    Returns [base_dir, base_dir_RSTRT1, ...] (base_dir always included).
    """
    base_dir = base_dir.rstrip('/')
    chain = [base_dir]
    n = 1
    while True:
        candidate = "{}_RSTRT{}".format(base_dir, n)
        if os.path.isdir(candidate):
            chain.append(candidate)
            n += 1
        else:
            break
    return chain


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

if len(sys.argv) < 2:
    print("usage: python3 plot_energies_combined.py <dir1> [<dir2> ...]")
    print("       (with a single <dir1>, _RSTRT1, _RSTRT2, ... are auto-discovered)")
    sys.exit(1)

if len(sys.argv) == 2:
    dirs = find_restart_chain(sys.argv[1])
    if len(dirs) < 2:
        print("no {}_RSTRT* restarts found; nothing to combine".format(sys.argv[1].rstrip('/')))
        sys.exit(1)
    print("auto-discovered restart chain: {}".format(dirs))
else:
    dirs = [d.rstrip('/') for d in sys.argv[1:]]

dirs = [d + '/' for d in dirs]

config = None
try:
    filename = "{}options.cfg".format(dirs[0])
    config = ConfigParser()
    config.read(str(filename))
except:
    print("failed to read config file. Please supply run suffix")
    raise

try:
    plot_energies(dirs, prefix='')
except Exception as e:
    print(e)
    print('failed to plot energies')
try:
    plot_zmodes(dirs)
except Exception as e:
    print(e)
    print('failed to plot zmodes')
try:
    plot_2energies(dirs, prefix='', config=config)
except Exception as e:
    print(e)
    print('failed to plot energies')
try:
    plot_invariant(dirs)
except Exception as e:
    print(e)
    print('failed to plot invariant')

try:
    plot_12_phase(dirs)
except Exception as e:
    print(e)
    print('failed to plot phase_12')

try:
    plot_zmodes(dirs, prefix='be')
except Exception as e:
    print(e)
    print('failed to plot zmodes')

try:
    plot_task(dirs, 'keff', logscale=False)
except Exception as e:
    print(e)
    print('failed to plot keff')

try:
    plot_task(dirs, 'ke_y', logscale=False)
except Exception as e:
    print(e)
    print('failed to plot ke_y')

try:
    plot_task(dirs, 'bhelicity', logscale=False)
except Exception as e:
    print(e)
    print('failed to plot bhelicity')

try:
    plot_task(dirs, 'be', logscale=True)
except Exception as e:
    print(e)
    print('failed to plot be')

try:
    plot_task(dirs, 'be_y', logscale=True)
except Exception as e:
    print(e)
    print('failed to plot be_y')

try:
    plot_task(dirs, 'ke_x', logscale=True)
except Exception as e:
    print(e)
    print('failed to plot ke_x')

try:
    plot_task(dirs, 'udiff', logscale=True)
except Exception as e:
    print(e)
    print('failed to plot udiff')

try:
    plot_task(dirs, 'Adiff', logscale=True)
except Exception as e:
    print(e)
    print('failed to plot Adiff')

try:
    plot_task(dirs, 'Adiff_hat', logscale=True)
except Exception as e:
    print(e)
    print('failed to plot Adiff_hat')

try:
    plot_task(dirs, 'omega', logscale=False)
except Exception as e:
    print(e)
    print('failed to plot omega')

try:
    plot_task(dirs, 'proj_A0')
except Exception as e:
    print(e)
    print('failed to plot projection A0')

try:
    plot_task(dirs, 'proj_norm_A0')
except Exception as e:
    print(e)
    print('failed to plot projection A0 normalized')

try:
    plot_task(dirs, 'be_y-be_z')
except Exception as e:
    print(e)
    print('failed to plot b-diff')

try:
    plot_task(dirs, 'Rossby')
except Exception as e:
    print(e)
    print('failed to plot Rossby')

try:
    plot_task(dirs, 'Rayleigh')
except Exception as e:
    print(e)
    print('failed to plot Rayleigh')

try:
    plot_task(dirs, 'AdotA_mean')
except Exception as e:
    print(e)
    print('failed to plot AdotA_mean')
