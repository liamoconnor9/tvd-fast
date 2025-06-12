# tvd-fast: Pseudospectral Simulations of the Taylor-vortex Dynamo

This repository contains code for performing kinematic and nonlinear dynamo simulations using the Dedalus pseudospectral framework. It accompanies the publication:

**"The Unsteady Taylor-Vortex Dynamo is Fast"**  
*by Liam O'Connor et al*  
*(Currently under review)*
## Overview
The code simulates magnetic field growth in an unsteady Taylor-vortex flow, demonstrating fast dynamo action (exponential field growth independent of magnetic diffusivity at low diffusivities)
- Distributed pseudospectral implementation using [Dedalus](http://dedalus-project.org/)
- Solves the induction equation for magnetic field evolution in a rotating-plane Couette flow configuration
## Usage
These accompanying scripts are compatible with supercomputers running the SLURM workload manager. Simulations can be performed as follows:
1. Customize the configuration ``options.cfg`` by specifying job/directory name, labelled "suffix"
2. Customize the configuration ``options.cfg`` by specifying the physical parameters (Re, Rm, resolution, timestep, domain size, etc.) and compute resources (# nodes, # cores, wall time)
3. Customize the configuration ``options.cfg`` by specifying the SOLVER and TEMPLATE:
	- set SOLVER='rpcf-mhd.py' and TEMPLATE='template.sh' to perform a single simulation. 
	- set SOLVER='ky_spectrum.py' and TEMPLATE='spectrum_template.sh' to solve over a spectrum of streamwise wavenumbers. Customize the range of wavenumbers by editing the for loop in spectrum_template.sh
4. On a SLURM system, deploy and queue the simulation by running ``bash slrun.sh -q`` to submit a job. The job name, job ID, and timestamp will be appended to a text file called ``pilot.txt``
		
	or...
		
5. In an interactive/development setting (SLURM or otherwise), a simulation can be deployed and run using ``bash slrun.sh``

