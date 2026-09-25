#!/usr/bin/env python3
from matplotlib.lines import Line2D
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

def load_xy(path):
    # Reads two-column CSV, ignoring comment lines starting with '#'
    data = np.loadtxt(path, delimiter=",", comments="#", dtype=float)
    # Ensure 2D even if there's only one row
    data = np.atleast_2d(data)
    # Sort by x in case the file isn't ordered
    idx = np.argsort(data[:, 0])
    return data[idx, 0], data[idx, 1]

alpha=0.8
ms=9
mew=0.5
series = [
    ("gr_0p18.csv", 0.18, dict(marker="s", linestyle="none", label=r"$k_y=0.18$", ms=ms, color="lime", alpha=alpha, mec="black", mew=mew)),
    ("gr_0p29.csv", 0.29, dict(marker="o", linestyle="none", label=r"$k_y=0.29$", ms=ms, color="magenta", alpha=alpha, mec="black", mew=mew)),
    ("gr_1p0.csv",  1.00, dict(marker="^", linestyle="none", label=r"$k_y=1.0$",  ms=ms, color="cyan", alpha=alpha, mec="black", mew=mew)),
]

plt.figure(figsize=(6.4, 4.0))

mss=8
for fname, ky, style in series:
    x, y = load_xy(fname)
    ogstyle = style.copy()
    if ky == 0.18:
        mask2 = x > 70
        Rm2 = x[mask2]
        gr2 = y[mask2]
        style['marker'] = 'o'
        style['ms'] = mss
        plt.plot(Rm2, gr2, **style)

        style = ogstyle
        style['ms'] = 10
        mask0 = x < 70
        Rm0 = x[mask0]
        gr0 = y[mask0]
        style['marker'] = '*'
        # style['alpha'] = 1
        plt.plot(Rm0, gr0, **style)
    elif ky == 0.29:
        style['marker'] = 'o'
        style['ms'] = mss
        plt.plot(x, y, **style)
    elif ky == 1.0:
        mask2 = x >= 1.5e4
        Rm2 = x[mask2]
        gr2 = y[mask2]
        style['marker'] = 'o'
        style['ms'] = mss
        plt.plot(Rm2, gr2, **style)
        style = ogstyle

        mask1 = (x < 1.5e4) & (x > 800)
        Rm1 = x[mask1]
        gr1 = y[mask1]
        style['marker'] = '^'
        plt.plot(Rm1, gr1, **style)
        style = ogstyle

        mask0 = x < 800
        Rm0 = x[mask0]
        gr0 = y[mask0]
        style['ms'] = 10
        style['marker'] = '*'
        # style['alpha'] = 1
        plt.plot(Rm0, gr0, **style)
    else:
        raise
        # plt.plot(x, y, **style)


# for fname, ky, style in series:
#     x, y = load_xy(fname)
#     plt.plot(x, y, **style)

# Axes styling to resemble the example
plt.xscale("log")
plt.xlabel(r"$\mathrm{Rm}$")
plt.ylabel("growth rate")

# Optional: tighten the view to a typical range like in your image
# (feel free to remove/adjust these if your data ranges differ)
# plt.ylim(0.0, 0.33)
plt.ylim(bottom=-0.04)

# Legend (place it so it doesn't cover data too much)
# plt.legend(frameon=False, loc="lower right")
color_handles = [
    Line2D([0], [0], marker='^', color='none', markerfacecolor='white', markersize=8, label=r"$n=1$"),
    Line2D([0], [0], marker='o', color='none', markerfacecolor='white', markersize=8, label=r"$n=2$"),
    Line2D([0], [0], marker='*', color='none', markerfacecolor='white', markersize=10, label=r"other"),
    Line2D([0], [0], marker='s', color='none', markerfacecolor='lime', markersize=8, label=r"$k_y=0.18$", alpha=alpha, mec="black", mew=0),
    Line2D([0], [0], marker='s', color='none', markerfacecolor='magenta', markersize=8, label=r"$k_y=0.29$", alpha=alpha, mec="black", mew=0),
    Line2D([0], [0], marker='s', color='none', markerfacecolor='cyan', markersize=8, label=r"$k_y=1.0$", alpha=alpha, mec="black", mew=0),
]

plt.legend(handles=color_handles, frameon=True, ncol=2, handletextpad=0.3)

plt.tight_layout()
figname = "gr_plot.png"
plt.savefig(figname, dpi=300)
print(figname)
figname = "gr_plot.pdf"
plt.savefig(figname, dpi=300)
print(figname)
# plt.show()
