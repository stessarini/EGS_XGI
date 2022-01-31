/*
###############################################################################
#
#   EGS_XGI Example_Usercode_Score_Energy header
#   An example usercode for using EGS_XGI features within EGSnrc with energy
#		scoring features (compare EGSnrc tutor7pp)
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

#ifndef _EXAMPLE_EXAMPLE_USERCODE_SCORE_ENERGY_H_
#define _EXAMPLE_EXAMPLE_USERCODE_SCORE_ENERGY_H_

#include <iostream>
#include <string.h>

#include "EGS_XGIApplication.h"
#include "egs_advanced_application.h"
#include "egs_interface2.h"

#include "egs_functions.h"	/**/
#include "egs_base_source.h"	/**/
#include "egs_timer.h"		/**/
#include "egs_rndm.h"
#include "egs_scoring.h"

class Example_Usercode_Score_Energy : public EGS_XGIApplication
{
private:

public:
	Example_Usercode_Score_Energy(int argc, char **argv);

	virtual ~Example_Usercode_Score_Energy();

	virtual int ausgab(int iarg);

	//virtual int run();

	//void ActivateAusgabCals();

	void startHistory(EGS_I64 this_case);

	void endHistory();

	void reportResults();

	//void finish();
protected:

	int m_nNumberOfRayTracingPropagations;
	int m_nNumberOfRayleighEvents;
//	int m_nStackSizeLastCall;
private:
	EGS_ScoringArray* m_pDepositedEnergy;
	int m_nNumberOfRegions;
	EGS_Float m_fScoredEnergyThisHistory;
	int m_nInteractionCheck;

};


#endif
