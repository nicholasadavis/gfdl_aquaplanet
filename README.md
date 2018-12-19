# gfdl_aquaplanet
Source code modifications for the GFDL Grey Radiation Aquaplanet Model

User Guide
============================

Users are referred to the Geophysical Fluid Dynamics Laboratory documentation for installing and running the Flexible Modeling System. 

Known issues
============================

Dec. 2018: In a pure sigma coordinate model, the stratospheric circulation in axisymmetric circulations with prescribed eddy forcings deviates substantially from the expected value. Running the model with hybrid-sigma coordinates, with the transition to pure pressure levels around 200 hPa, eliminates this behavior. A text file with namelist commands for a and b coefficients for a 40-layer hybrid-sigma model is included in this release. 

Dec. 2018: The shortwave ozone optical depth has a finite value at the model lid, resulting in extreme shortwave absorption at the model lid.

Installation
=======================================================================

I recommend you install a parallel compilation of the model if you wish to use this modified source code. Copy your source code, replace the relevant files found here, and compile the model. I have not verified that this modified code will function with any particular combination of compilers and on any particular architecture. I have verified it is stable running on 16 and 64 cores using the Intel 15.6.233, Intel MPI 5.0.3.048, NetCDF 4.4.0, and HDF5 Parallel 1.8.14 compilers on NOAA’s Theia architecture. Please contact me if you have questions about compilers or compiler options - I might be able to help.

Features
=======================================================================

Users are given the ability to prescribe either the zonal-mean eddy or mean flow. A description of the prescribed eddy technique is found in Davis and Birner (in review, you can contact me for the draft). The prescribed mean flow is a recent update and has not been rigorously tested. Regardless of which technique is used, this code is provided as-is with no guarantee of scientific performance. Users are also given the ability to enable the Levine and Schneider idealized ocean flux scheme, as well as an idealized shortwave absorption scheme meant to mimic ozone absorption and vertical diffusion (Davis and Birner, in review).

Users can prescribe the zonal-mean eddy tendencies from an eddy-permitting simulation to an axisymmetric simulation to isolate the response of the mean flow to forcings, or can prescribe the zonal-mean flow from an eddy-permitting simulation to another eddy-permitting simulation to isolate the response of eddies to forcings. Both of these will require the creation of an input data file in netcdf format. An example is included with this release. 

The simulation options provided here are very simple, and they can be used in ways not described here. 

Please reference Davis and Birner (in review) if you use this modified source code, and feel free to contact me with questions or ideas. 

Prescribing eddy tendencies to isolate the mean flow response to forcings
=======================================================================

Run an eddy-permitting simulation first to derive eddy tendencies for an axisymmetric simulation. The following variables must be saved: uu, dy(vv), dy(uv), dy(vq), dy(vt), dp(ov), dp(ou), dp(ot), dp(oq), ot, eddy_mass_div, eddy_mass_u, eddy_mass_v, eddy_u_stress, eddy_v_stress, eddy_t, eddy_q. I recommend also saving tdt_ls and qdt_ls (see Davis and Birner, in review). You only need to save these as average fields over the entire model run, there is no need to output to a higher time resolution. 
Combine these forcings into an input file in netcdf format. See the example input file included in this release. I have also included a Matlab script as an example post-processing file.
Run the axisymmetric simulation. The following spectral dynamics namelist parameters must be set: prescribe_fluxes = .true., input_file = /path/to/your/file.nc. You must also initialize the run with the spectral dynamics init namelist option do_zonal_perturbation = .false.

Prescribing the mean flow to isolate the eddy response to forcings
=======================================================================

Run an eddy-permitting simulation first. The following variables must be saved: u, v, t, q, ps. You only need to save these as average fields over the entire model run, there is no need to output to a higher time resolution.
Combine these into an input file. See the example input file included in this release. I have also included a Matlab script as an example post-processing file.
Run another eddy-permitting simulation. The following spectral dynamics namelist parameters must be set: prescribe_mean = .true., input_file = /path/to/your/file.nc. 

Enabling vertical diffusion
=======================================================================

Set do_diffusion = .true. in the spectral dynamics namelist; you can modify the diffusion coefficient via diffusion_coefficient, its default value is 5 (m^2/s).

Enabling shortwave psuedo-ozone absorption
=======================================================================

Set atm_abs = 0.008 and frac_ozone = 1.0 in the gray radiation namelist to replicate the ozone distribution in Davis and Birner (in review). If you are including tropospheric shortwave absorption, you will need to modify these values. E.g., for 10% tropospheric absorption, atm_abs = 0.108 and frac_ozone = .0741.

Other options
=======================================================================

The infrastructure for applying eddy tendencies does not assume the prescribed fields are zonally-averaged. It is therefore possible to use this infrastructure for forcings with zonal structure. 

Full description of changes
=======================================================================
The following files have been modified from their original state: spectral_dynamics.f90 spectral_init_cond.f90 spectral_initialize_fields.f90 simple_surface.f90 surface_flux.f90 gray_radiation.f90

spectral_dynamics.f90: modified to apply prescribed eddy tendencies to the model atmosphere to the u, v, T, ps, and q fields, or to apply prescribed zonal mean flow fields. Modified to save eddy flux convergences on-line in the diagnostics section. Modified to include a quadratic interpolation subroutine, fast eddy flux calculation subroutine, horizontal convergence subroutine, vertical convergence subroutine, and various other subroutines (zonal averaging, zonal mean flow forcing, etc.). The horizontal convergence subroutine employs the spectral method to calculate derivatives exactly, while the vertical convergence subroutine uses a centered difference method. At low vertical resolutions, this may be suboptimal. It is recommended that you run the model with at least 40 vertical levels so that the boundary layer, stratosphere, and vertical flux convergences are resolved. Includes the option to run with a constant free atmospheric diffusion for u, v, T, and q.

spectral_init_cond.f90: modified to accept namelist input to initialize the model as zonally-symmetric.

spectral_initialize_field.f90: modified to initialize the model as zonally-symmetric.

simple_surface.f90: modified to apply prescribed eddy contributions to surface drag and heat and momentum fluxes. Passes these fields to the surface_flux module. Includes modifications for prescribed SST input as well, though this is presented without description or validation. Use at your own risk. Modified to include the Levine and Schneider idealized meridional ocean heat flux scheme.

surface_flux.f90: modified to accept input eddy surface flux contributions and apply them to the surface transfer equations.

gray_radiation.f90: modified to create a shortwave absorption layer in the stratosphere that mimics the effects of ozone. This creates a stratospheric jet and, miraculously, separates the eddy-driven and subtropical jets in the troposphere. This minor fix for the climatology is recommended, but it has only been tuned to a surface longwave optical depth of 6 (2) at the equator (pole). I do not recommend running the model with fewer than 40 levels if you use this modification.
