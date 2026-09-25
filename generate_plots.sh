#!/bin/bash

cd ~/tvd-fast

suffix1="nl_ky1_Rm3e3_noiseIC"
suffix2="nl_ky1_Rm1p5e3_fromscratch"
suffix3="nl_Lz2pi_Rm1e3"

suffix4="nl_Lz2pi_Rm5e3_seed35"

# nl_Lz2pi_Rm5e3_seed35
python3 plot_oscillation_spectrum.py $suffix1
python3 plot_oscillation_spectrum.py $suffix2
python3 plot_oscillation_spectrum.py $suffix3
python3 plot_oscillation_spectrum.py floquet-Rm1e3
python3 plot_oscillation_spectrum.py floquet-Rm3e3

python3 plot_energies.py $suffix1
python3 plot_energies.py $suffix2
python3 plot_energies.py $suffix3
python3 plot_energies.py $suffix4
# python3 plot_energies.py nl_Lz2pi_Rm5e3_seed16

# python3 plot_oscillation_spectrum.py nl_Lz2pi_Rm1e3