/*
###############################################################################
#
#   EGS_XGI OpticalElement header
#   Base class for all optics components/elements, i.e. gratings, detectors,...
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
#ifndef _OPTICALELEMENT_H_
#define _OPTICALELEMENT_H_

#include "egs_vector.h"
#include "egs_input.h"

class EGS_Interpolator;
class EGS_BaseGeometry;
class OpticalElement
{
public:
	OpticalElement();
	virtual ~OpticalElement();
	virtual bool IsOnOpticalElement (EGS_Vector* i_pPosition);
	virtual bool ApplyOpticalElementRule();
	virtual void PrintSummary();
	virtual void ReportOpticalElement();
	//trigger precalculation of parameters that are constant for this shower call
	virtual void StartRayTracing(EGS_Float& i_fEnergy);
	//reset precalculated parameters
	virtual void EndRayTracing();
	virtual void ScoreSecondaryParticle();
	virtual int InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry);

	//Give the SplittingAlgorithm the permission to reset its pointer to the previous OpticalElement after all paths are transported
	virtual bool PermissionToResetPointer(int i_nStackSize);
	EGS_Float GetNumberOfTimesAppliedOpticalElemetRule();
	virtual void EndSimulation(int i_nNumberOfHistories);

	std::string GetName();

	double GetPhaseOfComplexNumber(EGS_Float i_fRealPart, EGS_Float i_fImaginaryPart);

protected:
	unsigned int m_nNumberOfTimesApplied;

	std::string m_sName;

	int m_nStackSizeBefore;
private:

};

#endif
