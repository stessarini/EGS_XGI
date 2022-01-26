/*
###############################################################################
#
#   EGS_XGI SplittingAlgorithm header
#   Interferometer base class, managing splitting of the optics components
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
#ifndef _SPLITTINGALGORITHM_H_
#define _SPLITTINGALGORITHM_H_
#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>

#include "OpticalElement.h"
#include "RefractiveIndexCalculator.h"

#include "egs_rndm.h"
#include "egs_interface2.h"
#include "egs_input.h"


class SplittingAlgorithm
{
public:

SplittingAlgorithm();

virtual ~SplittingAlgorithm();

SplittingAlgorithm(std::string i_sInputFile);

//Report relevant variables of the interferometer and OpticalElements
virtual void ReportSplittingAlgorithm();

//Report of OpticalElements after run
virtual void ReportSplittingSummary(int ncase);

//Do path splitting if necessary
virtual	void PotentialParticleSplitting();

//Reset Parameters, e.g. pointer to OpticalElement
virtual	void PotentialParameterReset();

int InitSplittingAlgorithm(EGS_Input* i_Input, std::string* i_sInputFileName, EGS_RandomGenerator* i_RandomNumberGenerator, EGS_BaseGeometry* i_pGeometry, EGS_Interpolator* i_gmfp, RefractiveIndexCalculator* i_pRefractiveIndexCalculator);

virtual int InitDerivedSplittingAlgorithm(EGS_Input* i_Input);

OpticalElement* GetPointerTo(std::string i_sName);

virtual void StartRayTracing(EGS_Float i_fInitialEnergy);

virtual void endHistory();

virtual void EndSimulation(int i_nNumberOfHistories);

virtual std::string WhichOpticalElement();

//return 0 if the splitting algorithm is in the correct state to start a new history
//otherwise erturn non-zero int
virtual int IsInInitialCondition();

protected:
vector<OpticalElement*> m_pSplittingObjects;
vector<OpticalElement*> m_pAllSplittingObjectsInInputfile; //All OpticalElements defined in the input file

};


#endif
