from configparser import ConfigParser
import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
import sys


def get_task(dir, task):
    files = glob.glob("{}scalars/*.h5".format(dir))
    files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    last_index = len(files)
    times_data = []
    data_data = []
    for file in files:
        with h5py.File(file, "r") as f:
            times_data += f['scales']['sim_time'][()].tolist()
            data_data += f['tasks'][task][()].ravel().tolist()

    return (times_data, data_data)


def get_mean_keff(window, parent):
    times, keff_data = get_task(parent + '/', 'keff')
    times = np.array(times)
    keff_data = np.array(keff_data)
    # window = 23.33
    t_max = times[-1]
    t_min = t_max - window
    mask = times >= t_min
    times_recent = times[mask]
    keff_recent = keff_data[mask]
    return np.mean(keff_recent)

suffixlist = "suffices_ky0p18.txt"
window = 23.33
outputfile = "keff_0p18.csv"

suffices = []
with open(suffixlist, "r") as f:
    for line in f:
        line = line.strip()  # remove newline and whitespace
        if not line.startswith("#") and line:  # skip comments and blank lines
            if not ' ' in line:
                suffices.append(line)
            else:
                suffices.append(line.split(' ')[0])

final = []
for suffix in suffices:
    try:
        filename = glob.glob("{}/options.cfg".format(suffix))[0]
        config = ConfigParser()
        config.read(str(filename))
        Rm = config.getfloat("parameters", "Rm")
        keff = get_mean_keff(window, suffix)
        final.append((Rm, keff**2 / Rm))
    except:
        print("failed to read config file. Bad suffix = {}".format(suffix))
        raise

with open(outputfile, "w") as f:
    for x, y in final:
        if x > 9e5:
            continue
        f.write(f"{x}, {y}\n")
        print(x)
        print(y)
print(outputfile)

# print(keff_recent[0])
# print(keff_recent[-1])

# print(times_recent[0])
# print(times_recent[-1])
# print(len(times_recent))