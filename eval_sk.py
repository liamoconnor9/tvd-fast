from configparser import ConfigParser
import numpy as np
# import dedalus.public as d3
import matplotlib.pyplot as plt
import h5py
import glob
import sys

def get_zmode_data(dir, prefix='ke'):

    files = glob.glob("{}/scalars/*.h5".format(dir))
    last_index = len(files)
    mode1 = []
    mode2 = []
    times = []
    for file in files:
        with h5py.File(file, "r") as f:
            times += f['scales']['sim_time'][()].tolist()
            mode1 += f['tasks']['{}_mode{}'.format(prefix, 1)][()].ravel().tolist()
            mode2 += f['tasks']['{}_mode{}'.format(prefix, 2)][()].ravel().tolist()
    return (times, mode1, mode2)

def get_kill_time(dir):
    times, mode1, mode2 = get_zmode_data(dir)
    zipped_lists = zip(times, mode1, mode2)

    sorted_zipped_lists = sorted(zipped_lists)
    new_times, new_mode1, new_mode2 = zip(*sorted_zipped_lists)

    times = list(new_times)
    mode1 = list(new_mode1)
    mode2 = list(new_mode2)

    time_np = np.array(times)
    mode1_np = np.array(mode1)
    
    if len(mode1_np) == 0:
        return None
    
    final_value = mode1_np[-1]
    
    abs_diff = np.abs(mode1_np - final_value)
    tolerance = 1e-5
    not_constant_indices = np.where(abs_diff > tolerance)[0]
    
    if not_constant_indices.size == 0:
        return time_np[0]
    
    constant_start_index = not_constant_indices[-1] + 1
    
    if constant_start_index >= len(time_np):
        if times[-1] < 9999:
            print(times[-1])
            return -times[-1]
        else:
            return None
    else:
        return time_np[constant_start_index]

    # end_time = np.max(times)
    # def check_window(win_arg):
    #     mode1_final_val = mode1[-1]
    #     filter_bool = times > end_time - win_arg
    #     end_mode1 = np.array(mode1)[filter_bool]
    #     return np.allclose(end_mode1 - mode1_final_val, 0)

    # window = 5
    # killed = check_window(window)
    # if not killed:
    #     return None
    # else:
    #     for i in range(window, int(end_time)):
    #         killed_i = check_window(i)
    #         if not killed_i:
    #             return round(end_time - i)

suffixlist = "sk-stats.txt"
suffices = []
with open(suffixlist, "r") as f:
    for linenum, line in enumerate(f):
        line = line.strip()  # remove newline and whitespace
        if not line.startswith("#") and line:  # skip comments and blank lines
            if not ' ' in line:
                suffices.append(line)
            else:
                suffices.append(line.split(' ')[0][:-1])

final = []
for suffix in suffices:
    try:
        filename = glob.glob("{}/options.cfg".format(suffix))[0]
        config = ConfigParser()
        config.read(str(filename))
        # print(suffix)
        kt = get_kill_time(suffix)
        if kt > 9995:
            kt = None
        elif kt < 0:
            kt = 'IDK'
        final.append((suffix, kt))
    except:
        print("failed to read config file. Bad suffix = {}".format(suffix))
        raise

outputfile = "kill_times.csv"
with open(outputfile, "w") as f:
    for x, y in final:
        f.write(f"{x}, {y}\n")
print(outputfile)