/*
###############################################################################
#
#   EGS_XGI global variables, additional photon stack variables
#   Class to handle the decrement \delta of the real part of the refractive
#   index: n = 1 - \delta
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
*/

#include "xgi_global_variables.h"

EGS_Float phase_array[MXSTACK];
EGS_Float* e_fPhase = phase_array;

EGS_Float norm_array[MXSTACK];
EGS_Float* e_fLogNorm = norm_array;

bool primary_array[MXSTACK] = {false};
bool* e_bPrimary = primary_array;

const EGS_Float ec_fPi = 3.141592653589793238463;

const EGS_Float ec_fEnergyToWaveLength = 1.23984187541999*1e-10; /*\lambda[cm] = c_fEnergyToWaveLength/E[MeV]*/

const std::string ec_slibVersion = "EGS_XGI_v1.0";
