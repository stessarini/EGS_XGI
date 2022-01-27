#!/bin/bash
###############################################################################
#
#   Bash script to run all sets of Talbot carpet simulations for the comparison
#   between uniform (Huygens) and Fourier splitting.
#   Copyright (C) 2020  ETH Zürich
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
#   Author:     Stefan Tessarini
#
#
#
###############################################################################
#set the following variables:
export EGS_HOME=/path/to/egs_xgi_home/
#if not using default:
export HEN_HOUSE=/path/to//HEN_HOUSE/
export EGS_CONFIG=/path/to/config/file
export my_machine=my_machine_name

echo "#####################################"
echo "run__TC_Fourier_1e6"
./run__TC_Fourier_1e6.sh
echo "#####################################"
echo "run__TC_Fourier_4e6"
./run__TC_Fourier_4e6.sh
echo "#####################################"
echo "run__TC_Fourier_8e6"
./run__TC_Fourier_8e6.sh
echo "#####################################"
echo "run__TC_Huygens_1e6"
./run__TC_Huygens_1e6.sh
echo "#####################################"
echo "run__TC_Huygens_5e5"
./run__TC_Huygens_5e5.sh
echo "#####################################"
pwd
