#!/bin/bash
N=$1
sed -i "s/^suffix=.*/suffix='skd${N}_Ro3p5_Rm1p2e3_Ny16_Lz2pi'/" options.cfg
sed -i "s/^seed=.*/seed=${N}/" options.cfg
bash slrun.sh -q
