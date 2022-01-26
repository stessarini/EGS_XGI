/*
###############################################################################
#
#   EGS_XGI DefaultSplittingAlgorithm header
#   Default interferometer class, managing splitting of the optics components
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

#ifndef _DEFAULTSPLITTINGALGORITHM_H_
#define _DEFAULTSPLITTINGALGORITHM_H_
#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>

#include "SplittingAlgorithm.h"

#include "egs_input.h"

#include "egs_rndm.h"



class DefaultSplittingAlgorithm : public SplittingAlgorithm
{
public:
	DefaultSplittingAlgorithm();

	~DefaultSplittingAlgorithm();

	void ReportSplittingAlgorithm();

	void ReportSplittingSummary(int ncase);

	void PotentialParticleSplitting();

	void PotentialParameterReset();

	virtual void endHistory();

	virtual void EndSimulation(int i_nNumberOfHistories);

	virtual int InitDerivedSplittingAlgorithm(EGS_Input* i_Input);

	virtual std::string WhichOpticalElement();

	virtual int IsInInitialCondition();

protected:
	//The last OpticalElement used for path splitting:
	vector<OpticalElement*>::iterator m_pLatestActiveSplittingObject;

	bool m_bCoherentSource;

private:

};

#endif
