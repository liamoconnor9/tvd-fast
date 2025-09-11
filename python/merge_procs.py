import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from dedalus.extras import plot_tools
from dedalus.tools import post
import sys


post.merge_analysis('Ly22_kin_Rm22p0_RSTRT2/checkpoint', cleanup=True)
# post.merge_process_files('dataSlices/restart-1Re150.0__Ro3.0/slices', cleanup=True)
