#!/usr/bin/env python3
from matplotlib.lines import Line2D
import numpy as np
import matplotlib.pyplot as plt
import os
path = os.path.dirname(os.path.abspath(__file__))
import matplotlib.pyplot as plt
# plt.rcParams.update({'font.size': 14})
import sys
# plt.style.use('dark_background')

def load_xy(file_name):
    data = []
    with open(file_name) as f:
        for line in f:
            a, b = line.strip().split(",")
            a, b = a.replace(" ", ""), b.replace(" ", "")
            data.append((a, b))    
    col1, col2 = zip(*data)
    col1, col2 = list(col1), list(col2)
    return (col1, col2)

plt.figure(figsize=(4, 3))
datax, datay = load_xy(path + '/kill_times.csv')

data_dict = {}
for i in range(len(datax)):
    name = datax[i]
    val = datay[i]
    if not "RSTRT" in name:
        data_dict[name] = eval(val)
        if data_dict[name] == None:
            data_dict[name] = 10000
    else:
        spl_name = name.split("_RSTRT")
        og_name = spl_name[0]
        value = eval(val)
        if value == None:
            data_dict[og_name] += 10000
        else:
            data_dict[og_name] += value
        # if val == None:
        # else:
        #     data_dict[name] = float(val)
lifetimes = list(data_dict.values())
print(lifetimes)
T_bar = np.mean(lifetimes)
# std = np.var(lifetimes)
print('T_bar = {}'.format(T_bar))
T_max = 3.6e4
T = np.linspace(0, T_max, 10000)
p_T = 1/T_bar * np.exp(-T / T_bar)

bins=6
bins = [0, 6e3, 1.2e4, 1.8e4, 2.4e4, 3e4, 3.6e4]
plt.hist(lifetimes, density=True, bins=bins, color='magenta', label='simulations')
plt.plot(T, p_T, color='deepskyblue', linewidth=5, label='best fit', alpha=1, linestyle='dotted')
plt.xlabel(r'$t_{\rm{life}}$')
plt.ylabel(r'$p(t_{\rm{life}})$')
plt.title(r"$\rm{Rm}=5000$")
plt.legend(frameon=False)
plt.tight_layout()
figname = path + "/lifetimes.png"
plt.savefig(figname, dpi=300)
print(figname)
figname = path + "/lifetimes.pdf"
plt.savefig(figname, dpi=300)
print(figname)
# plt.show()
