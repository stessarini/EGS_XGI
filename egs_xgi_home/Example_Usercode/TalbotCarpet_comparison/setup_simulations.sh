#!/bin/bash

###############################################################################
#
#   Bash script to set up the Talbot carpet simulations for the comparison
#   between uniform (Huygens) splitting and Fourier splitting
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


#Make sure a suitable python version is installed
#Dependent on your python version/system setup change python to python3,
# python3.8,.. if needed

mkdir "Fourier_1e6_"{5..43..2}
python create_input_files_TalbotCarpet_Fourier_1e6.py

mkdir "Fourier_4e6_"{5..43..2}
python create_input_files_TalbotCarpet_Fourier_4e6.py

mkdir "Fourier_8e6_"{5..43..2}
python create_input_files_TalbotCarpet_Fourier_8e6.py

mkdir "Huygens_1e6_"{0..6}
python create_input_files_TalbotCarpet_Huygens_1e6.py

mkdir "Huygens_5e5_"{0..6}
python create_input_files_TalbotCarpet_Huygens_5e5.py
