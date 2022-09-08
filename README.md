# 3d_mix_lck
**WARNING: CODE NOT FUNCTIONING, STILL IN DEVELOPMENT**
## Overview 

**Solver of the density-locked extend Gross-Pitaevskii (GP) equation in three-dimensions (3D)**

The code time steps using the Split-Step Fourier Method (SSFM), with both imaginary and real time. It relies on several libraries and packages:
* `FFTW3` - A standard FFT library
* `OpenMP` - Shared memory parallelisation library
* `HDF5` - Data format, used for outputting data here
* `json-fortran` - A JSON API for Fortran, used for input files for each simulation here

The code is structured as follows:

### `gp_lck.f03`
The `gp_lck` program is the main program script of the code which loads in input data (see `config.json`) and calls other functions and subroutines for setting up grid, time-stepping, etc.

### `FFTW.f03`
The `FFTW` module is simply to load FFTW into the code.

### `grid.f03`
The `grid` module contains three functions for setting up the grids and differential operators:
1) `space_grid(Nr,dr)` - The spatial 1D grid arrays (`Nr` - the number of grid points, `dr` - the spatial step size)
2) `mom_grid(Nr,dr)` - The momentum space 1D grid arrays 
3) `exp_lap(kx,ky,kz,dt)` -The Laplacian different operator written in the form `exp(-0.5*dt*(kx^2 + ky^2 + kz^2))`, as this is the form that is called within the SSFM time-stepping, so by defining this operator initially, there is a slight reduction in operations (`kx`, `ky` and `kz` are the three momentum space arrays, `dt` the time-step)

### `init.f03`
The `init` module (currently) contains one function for initialising the imaginary time-stepping (there is an aim to extend this to read in a `.h5` file as the initial wavefunction for real time):
1) `init_wa(x,y,z,init_type,gauss_sig)` - The initial wavefunction form (`x`, `y` and `z` are the three real space arrays, `init_type` allows for the selection of: a Gaussian initial density profile (`init_type = 1`); or a Super_Gaussian initial density profile (`init_type = 2`))

### `time.f03`
The `time` module contains two subroutines and one function for time-stepping the right-hand-side of the GP equation:
1) `ssfm(psi,dk2,t_steps,t_save,dt,dx,dy,dz,Nlck,mu,im_real)` - The SSFM is employed here and can be used for both imaginary and real time (though the time-step `dt` must be defined to be real or imaginary before being passed to this subroutine). The main structure of this method is:
* Half potential step in real space
* Transform to Fourier space and compute the kinetic energy term
* Transform back to real space and compute the other half potential step
* Finally, additional calculations include chemical potential and renormalisation if in imaginary time
This subroutine is also responsible for the in-time data outputting (to `.h5` format), and for creating the FFTW plans (`psi` is the wavefunction, `dk2` is the Laplacian operator in exponential form, `t_steps` number of time-steps, `t_save` number of steps between writing data, `Nlck` is the effective atom number for the density-locked system, `mu` is the chemical potential, and `im_real` denotes imaginary (`im_real = 0`) or real (`im_real = 1`) time)
2) `renorm(psi,dx,dy,dz,Nlck)` - This subroutine renormalises the wavefunction (computed at each time-step in imaginary time)
3) `chem_pot(psi,psi_k,dk2,plan_back,Nx,Ny,Nz,dt)` - The chemical potential is calculated at each imaginary time-step, to aid in imaginary time convergence (`psi_k` is the Fourier transform of the wavefunction, and `plan_back` is the Inverse Fourier Transform plan)

### `rhs.f03`
The `rhs` module contains two subroutines:
1) `V_rhs(psi,mu,dt,Nx,Ny,Nz)` - The potential half-step.
2) `T_rhs(psi_k,dk2,Nx,Ny,Nz)` - The kinetic energy step in momentum space.

## Compiling and Running
Once the requisite libraries are loaded and located on the system, the next step is to setup up the `config.json` (the naming convention the code expects) file from the `default.config.json` file. 

With the codes ready, the paths to the necessary libraries should be added to the associated Makefile (`makefile_pc` for local machines and `makefile_hpc` for remote clusters). Then copy the desired Makefile to `makefile` and run `make` to compile. **NOTE!** These Makefiles are highly simplified and will be formalised soon. If running locally then the code can be ran with the `./gplck` command, whilst for running on a cluster, an example SLURM job script is given in `3d_drop.sh`.

## Future...
Extension work will be undertaken on this code to include:
1) Processing scripts
2) Experimental to theoretical scripts to setup the `config.json` input file from experimental parameters
3) Ability to read in ground state wavefunction from `*.h5` file to propagate in real time

