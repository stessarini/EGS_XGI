# EGS_XGI
EGS_XGI is a c++ extension for EGSnrc (doi:10.4224/40001303) for Monte Carlo (MC) simulation of X-ray grating interferometry as published in S. Tessarini et al. . 
EGS_XGI adds diffraction effects to EGSnrc MC photon transport enabling simultaneous simulation of interference and scattering effects such as deposited energy. If properly installed EGS_XGI does not replace nor influence EGSnrc installations similar to a c++ EGSnrc usercode. For performing conventional EGSnrc simulations after the installation of EGS_XGI simply use the EGSnrc Makefiles and classes delivered by EGSnrc.



## License:
EGS_XGI is distributed as free software under the terms of the GNU Affero General Public License. Please consider the license before downloading the software.


## Prerequisites:
EGS_XGI can be installed as plug-in/extension for EGSnrc on Linux (WIN and MAC currently not supported) systems. To use EGS_XGI an EGSnrc installation is required (see: doi:10.4224/40001303), which has the minimal requirements:
 - Fortran compiler 
 - C compiler 
 - C++ compiler 
 - the GNU make utility


## Installation:
After successful installation of EGSnrc. Download the code or use git to clone the repository.
Unpack the code if necessary.
The extension comes with two main folders:
- EGS_XGI contains all extensions for EGSnrc (such as a base class for EGS_XGI applications, available optics components, makefiles, etc.)
- egs_xgi_home (the equivalent to EGS_HOME): a space for EGS_XGI usercodes (equivalent to EGSnrc usercodes). It contains:
  - bin: folder for the executable
  - Example_Usercode: an example usercode – a simple EGS_XGI application
  - Example_Usercode_Score_Energy: an EGS_XGI usercode that scores deposited energy on region basis similar to the EGSnrc tutor2pp usercode.

For installation a additional folder with the name $my_machine has to be created inside egs_xgi_home/bin/.
The usercodes can be compiled by:
- cd into the corresponding usercode folder, e.g., egs_xgi_home/Example_Usercode/
- set the following environment variables for compilation and simulation:
  EGS_HOME, HEN_HOUSE, EGS_CONFIG, XGI_PATH (only during compilation), and my_machine
  To do that open the Makefile set the environment variables (if they are not set globally): 
	  - EGS_HOME has to point to the egs_xgi_home folder
	  - HEN_HOUSE, my_machine, and EGS_CONFIG according to the EGSnrc installation that you are using
	  - XGI_PATH has to be set to EGS_XGI/
- then use the make command in the terminal.

## Run simulations:
A few example input files to reproduce the results in S. Tessarini et al. Are provided in the subfolders in Example_Usercode/ and  Example_Usercode_Score_Energy/. The simulations can be executed in the terminal by e.g. cd into your egs_xgi_home:
cd your_path_to/egs_xgi_home/Example_Usercode/DoubleSlit/

open the run script, e.g. with vi:
vi run__DoubleSlit

change the environment variables as done for the Makefile (EGS_HOME, HEN_HOUSE, EGS_CONFIG, and my_machine, except for XGI_PATH, which is not needed during the MC run). After successful compilation (see above) the simulation can be executed with the run script:
./run__DoubleSlit

## Example inputs:
The minimum input to generate the results in S.Tessarini et al. . Is included in the form of input files  where feasible. Cases that require more input files such as Talbot carpets are set up by python scripts that generate the input files. When using setup scripts such as setup_simulations.sh in egs_xgi_home//Example_Usercode/TalbotCarpet_comparison/ make sure that the python executable is defined in your PATH environment variable. Otherwise change the setup script accordingly.

## Python scripts:
The python scripts were executed with Python 3.8.10 and require the numpy and matplotlib.pyplot packages.


## Issues:

## Contributing:
