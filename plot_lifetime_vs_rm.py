#!/usr/bin/env python3
import glob
import os
import re
from configparser import ConfigParser

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

path = os.path.dirname(os.path.abspath(__file__))

# Each pattern is glob'd (relative to this file) to find every realization
# (seed) at a given Rm. Rm itself is read from each realization's
# options.cfg rather than parsed out of the directory name. Add new patterns
# here as more data comes in, rather than copying the whole script.
PATTERNS = [
    'skd*_Ro3p5_Rm8e2_Ny16_Lz2pi',
    'skd*_Ro3p5_Rm1e3_Ny16_Lz2pi',
]


def get_param(sim_dir, name):
    cfg = ConfigParser()
    cfg.read(os.path.join(sim_dir, 'options.cfg'))
    return cfg.getfloat('parameters', name)


def lifetime(sim_dir):
    """Last sim_time recorded in this run's scalar output, i.e. how long it
    ran before the dynamo died or it hit stop_sim_time. rpcf-mhd.py exits
    the instant it detects either condition, so the last scalar sample is
    the run's lifetime (to within one scalars_sim_dt)."""
    files = glob.glob(os.path.join(sim_dir, 'scalars', '*.h5'))
    t_max = None
    for f in files:
        with h5py.File(f, 'r') as h:
            t = h['scales']['sim_time'][()].ravel()
        if t.size:
            t_max = t.max() if t_max is None else max(t_max, t.max())
    return t_max


def base_name(sim_dir):
    # collapse "<base>_RSTRT<n>" restarts back onto "<base>" so a
    # restarted run's lifetime is the sum of its segments
    return re.sub(r'_RSTRT\d+$', '', os.path.basename(sim_dir.rstrip('/')))


def collect_lifetimes(pattern):
    dirs = sorted(glob.glob(os.path.join(path, pattern)))
    if not dirs:
        raise ValueError("no directories matched pattern '{}'".format(pattern))

    segment_sum = {}   # base run name -> summed lifetime across restarts
    any_censored = {}  # base run name -> True if any segment never died
    rm_values = set()
    for d in dirs:
        rm_values.add(get_param(d, 'Rm'))
        t_stop = get_param(d, 'stop_sim_time')
        t_life = lifetime(d)
        hit_stop = t_life is None or t_life >= t_stop - 1.0
        contribution = t_stop if t_life is None else t_life

        name = base_name(d)
        segment_sum[name] = segment_sum.get(name, 0.0) + contribution
        any_censored[name] = any_censored.get(name, False) or hit_stop

    if len(rm_values) > 1:
        raise ValueError("pattern '{}' matches runs with different Rm: {}".format(
            pattern, sorted(rm_values)))
    rm = rm_values.pop()

    lifetimes = np.array(list(segment_sum.values()))
    n_censored = sum(any_censored.values())
    if n_censored:
        print("warning: {} of {} runs in '{}' never died (hit stop_sim_time); "
              "their lifetime is a lower bound".format(n_censored, len(lifetimes), pattern))
    return rm, lifetimes


def main():
    results = sorted((collect_lifetimes(p) for p in PATTERNS), key=lambda r: r[0])

    rms = np.array([rm for rm, _ in results])
    means = np.array([lt.mean() for _, lt in results])
    stds = np.array([lt.std() for _, lt in results])

    for rm, (_, lt) in zip(rms, results):
        print("Rm={:g}: n={}, mean={:.1f}, std={:.1f}".format(rm, len(lt), lt.mean(), lt.std()))

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.errorbar(rms, means, yerr=stds, fmt='o-', color='darkviolet',
                capsize=4, markersize=6, label=r'$\langle t_{\rm{life}} \rangle \pm \sigma$')
    ax.set_xlabel(r'$\rm{Rm}$')
    ax.set_ylabel(r'$t_{\rm{life}}$')
    ax.set_ylim((0, None))
    ax.legend(frameon=False)
    plt.tight_layout()

    output_dir = os.path.join(path, 'lifetimes')
    os.makedirs(output_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        figname = os.path.join(output_dir, 'lifetime_vs_rm.{}'.format(ext))
        plt.savefig(figname, dpi=300)
        print(figname)
    plt.close(fig)


if __name__ == '__main__':
    main()
