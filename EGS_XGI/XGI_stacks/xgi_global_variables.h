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

#ifndef XGI_STACK_H_
#define XGI_STACK_H_

#include <array_sizes.h>
#include "egs_config1.h"
#include <string>


extern EGS_Float* e_fPhase;

extern EGS_Float* e_fLogNorm;

extern bool* e_bPrimary;

extern const EGS_Float ec_fPi;

extern const EGS_Float ec_fEnergyToWaveLength; /*\lambda[cm] = c_fEnergyToWaveLength/E[MeV]*/

extern const std::string ec_slibVersion;

#endif
