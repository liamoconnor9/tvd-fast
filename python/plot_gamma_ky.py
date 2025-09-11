import os
path = os.path.dirname(os.path.abspath(__file__))
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
data = []
import sys

# python3 plot_gamma_ky.py Rm36p86_3 Rm150 Rm861
# Rms = [36.9, 150, 861]

Rms = [36.9, 150, 861, 3000, 10000]
Rms_labels = ["36.9", "150", "861", "3000", r"$10^4$"]
# python3 plot_gamma_ky.py rm36p9_kin_spectrum_T466p6 rm150_kin_spectrum_T466p6 rm861_kin_spectrum_T466p6 rm3000_kin_spectrum_T466p6 rm1e4_kin_spectrum_T466p6
suffices = sys.argv[1:]
colors = ['lime',
  'cyan',
  'magenta',
  'blue',
  'orangered',
  'purple']
golden_mean = (np.sqrt(5)-1.0)/2.0
ratio = 1.0
# ratio = 1.075
fig_width = 5 # in column width
# fig_width = 7.1 # in page width
fig_height = golden_mean * fig_width
fig_size =  [fig_width, fig_height]
# plt.figure(figsize=fig_size)
plt.rcParams.update({'font.size': 16})
fig = plt.figure(figsize=(6, 4))

for ind, suffix in enumerate(suffices):
    ky_list = []
    gr_list = []

    isFirstLine = True
    with open('{}/{}/root_output0.txt'.format(path, suffix), 'r') as file:
        while True:
            temp = file.readline().replace('\n', '')
            if isFirstLine:
                isFirstLine = False
                continue
            if len(temp) == 0:
                break
            nums_strs = temp.split(', ')
            ky_list.append(float(nums_strs[0]))
            gr_list.append(float(nums_strs[1]) / 2)
            # print(temp)


    plt.plot(ky_list, gr_list, color=colors[ind], label = 'Rm={}'.format(Rms_labels[ind]), linewidth=2.5)
plt.hlines([0.0], xmin=0, xmax=1, color="grey", linestyle='dashed', linewidth=2.5)
plt.xlabel(r"$k_y$")
plt.ylabel('growth rate')
plt.xlim(0, 1)
# fig.subplots_adjust(top=0.9)
plt.legend(loc=(-0.15, 1.0), ncol=3, handlelength=1, frameon=False, fontsize=15)
# plt.legend(loc='upper center', bbox_to_anchor=(0.45, 1.35), ncol=3, handlelength=1, frameon=False)

# plt.legend()
# plt.title('Rm=36.86')
plt.tight_layout()
figname = "{}/{}/gamma_v_ky.png".format(path, suffix)
plt.savefig(figname)
print(figname)
figname = "{}/{}/gamma_v_ky.pdf".format(path, suffix)
plt.savefig(figname)
print(figname)