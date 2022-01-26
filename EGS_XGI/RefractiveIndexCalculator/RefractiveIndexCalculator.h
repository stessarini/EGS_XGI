/*
###############################################################################
#
#   EGS_XGI RefractiveIndexCalculator header
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

#ifndef _REFRACTIVEINDEXCLCULATOR_H_
#define _REFRACTIVEINDEXCLCULATOR_H_
#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>

#include "egs_base_geometry.h"

class RefractiveIndexCalculator
{
public:
	RefractiveIndexCalculator();
	~RefractiveIndexCalculator();
	void ImportRefractiveIndex(string i_sRefFile, EGS_BaseGeometry* i_pSimulationGeometry);
	EGS_Float GetRefractiveIndexDelta(int i_nMediumIndex, EGS_Float i_fEnergy);
	EGS_Float GetRefractiveIndexDeltaFromWavelength(int i_nMediumIndex, EGS_Float i_fWavelength);


private:
	//The refractive index data used to estimate the refractive index. Assume $\lambda^2$-dependence
	std::vector<EGS_Float> m_fRefractiveIndexData;
};


#endif
