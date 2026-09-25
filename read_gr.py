from scipy import stats
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


def get_gr(parent):
    times, bey_data = get_task(parent + '/', 'be_y')
    times = np.array(times)
    bey_data = np.array(bey_data)
    # window = 23.33
    t_max = times[-1]
    window = 0.9*t_max
    t_min = t_max - window
    mask = times >= t_min
    times_recent = times[mask]
    bey_recent = bey_data[mask]

    x, y = times_recent, np.log(bey_recent)
    slope, intercept, r, p, std_err = stats.linregress(x, y)
    print("parent={}; rsqrd={}".format(parent, r**2))
    return slope / 2

suffixlist = "suffices_ky0p18.txt"
outputfile = "gr_0p18.csv"

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
        gr = get_gr(suffix)
        final.append((Rm, gr))
    except:
        print("failed to read config file. Bad suffix = {}".format(suffix))
        raise

with open(outputfile, "w") as f:
    for x, y in final:
        f.write(f"{x}, {y}\n")
        # print(x)
        # print(y)
print(outputfile)