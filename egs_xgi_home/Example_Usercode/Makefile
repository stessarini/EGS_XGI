###############################################################################
#
#   Default EGS_XGI C++ makefile for the example usercode. Derived from EGSnrc
#   tutor7pp application
#
#
#   This file is part of the EGS_XGI - an X-ray grating interferometry
#   extension for EGSnrc.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as published
#   by the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
###############################################################################
#
#   Changes made by:     Stefan Tessarini
#
#
#
###############################################################################

# Original header:
###############################################################################
#
#  EGSnrc makefile to build the tutor7pp application
#  Copyright (C) 2015 National Research Council Canada
#
#  This file is part of EGSnrc.
#
#  EGSnrc is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Affero General Public License as published by the
#  Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  EGSnrc is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
#  more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with EGSnrc. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
#
#  Author:          Iwan Kawrakow, 2005
#
#  Contributors:
#
###############################################################################
#
#  As long as your own application is derived from EGS_AdvancedApplication and
#  there is a single C++ source file named $(USER_CODE).cpp, you can simply
#  use this Makefile by copying it into the application directory and changing
#  the USER_CODE variable below to the name of your application.
#
###############################################################################

#The following environment variables have to be defined:
#EGSnrc specific:
#  Only needed if they are not already set globally
#     HEN_HOUSE
#     EGS_CONFIG
#     my_machine
#
# EGS_HOME has to be temporarily reset for EGS_XGI simulations. (Both
# compilation and run-time)
# EGS_HOME has to contain 2 folders <bin> and one with the usercode name
# bin/ has to contain a folder with the name of $my_machine
# Extension library specific:
# XGI_PATH has to point to the location of the downloaded EGS_XGI folder

EGS_HOME := /path/to/egs_xgi_home/
HEN_HOUSE := /path/to/EGSnrc/HEN_HOUSE/
#usually the EGS_CONFIG file is located at: $(HEN_HOUSE)/specs/$(my_machine).conf
EGS_CONFIG := /EGSnrc/config/file
#sometimes required:
#ABS_EGSPP := /$(HEN_HOUSE)/egs++/
XGI_PATH := /path/to/EGS_XGI/
my_machine := gcc750

include $(EGS_CONFIG)
include $(SPEC_DIR)egspp1.spec
include $(SPEC_DIR)egspp_$(my_machine).conf


opt += -I$(XGI_PATH)/XGI_stacks/
opt += -I$(XGI_PATH)/EGS_XGIApplication/

opt += -I$(XGI_PATH)/Diffraction/OpticalElement/

opt += -I$(XGI_PATH)/Diffraction/OpticalElement/DetectorObject/DetectorZPlane/
opt += -I$(XGI_PATH)/Diffraction/OpticalElement/DetectorObject/DetectorZPlaneMem/

opt += -I$(XGI_PATH)/Diffraction/OpticalElement/SplittingObject/
opt += -I$(XGI_PATH)/Diffraction/OpticalElement/SplittingObject/FourierSeriesZPlane/
opt += -I$(XGI_PATH)/Diffraction/OpticalElement/SplittingObject/BinarySourceGratingFourierZPlane/
opt += -I$(XGI_PATH)/Diffraction/OpticalElement/SplittingObject/HuygensZPlane/


opt += -I$(XGI_PATH)/Diffraction/SplittingAlgorithm/
opt += -I$(XGI_PATH)/Diffraction/SplittingAlgorithm/DefaultSplittingAlgorithm/

opt += -I$(XGI_PATH)/RefractiveIndexCalculator/



opt += -std=c++0x

# Specify the name of the user code.
# The name of the executable is determined from this variable.
#
USER_CODE = Example_Usercode

# The following can be used to add user macros and mortran subroutines.
# The file(s) specified here are added after egsnrc.macros, machine.macros
# and egs_c_interface2.macros but before any files that have
# executable code.
#
EGSPP_USER_MACROS =

extra_user_files := xgi_global_variables

# Specify from which base class this application is being derived.
# This has the effect of automatically compiling the base application
# class and including it into the list of object files.
#
EGS_BASE_APPLICATION = egs_advanced_application

# Specify the set of mortran sources used.
# Here we simply use the standard set defined in $HEN_HOUSE/specs/egspp1.spec
#
CPP_SOURCES = $(C_ADVANCED_SOURCES)

# Specify here other header files that your user code depends upon.
#
other_dep_user_code =

include $(XGI_PATH)/Makefiles/cpp_makefile_xgi_application

egs_config1_path = $(HEN_HOUSE)/lib/$(my_machine)/

xgi_global_variables_$(my_machine).$(obje): $(XGI_PATH)/XGI_stacks/xgi_global_variables.cpp $(XGI_PATH)/XGI_stacks/xgi_global_variables.h $(egs_config1_path)/egs_config1.h
