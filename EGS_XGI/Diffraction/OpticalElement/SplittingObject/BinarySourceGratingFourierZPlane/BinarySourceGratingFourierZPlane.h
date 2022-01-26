/*
###############################################################################
#
#   EGS_XGI BinarySourceGratingFourierZPlane header
#	  Binary source grating orthogonal to the z-axis.
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
#ifndef _BINARYSOURCEGRATINGFOURIERZPLANE_H_
#define _BINARYSOURCEGRATINGFOURIERZPLANE_H_

#include <fstream>
#include <vector>

#include "egs_rndm.h"

#include "SplittingObject.h"

class BinarySourceGratingFourierZPlane : public SplittingObject
{
public:

	BinarySourceGratingFourierZPlane(std::string i_sName, EGS_RandomGenerator* i_pRandomNumberGenerator);
	BinarySourceGratingFourierZPlane(std::string i_sName);
	~BinarySourceGratingFourierZPlane();

	virtual bool IsOnOpticalElement (EGS_Vector* i_pPosition);
	virtual bool ApplyOpticalElementRule();
	virtual void ReportOpticalElement();
	virtual int InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry);
	virtual void PrintSummary();


private:
	EGS_Float m_fPosition;

	EGS_Float m_fHoleSizeInCM;
	EGS_Float m_fPeriodicityInCM;

	EGS_Float m_fTransmissionNorm;

	EGS_Float m_fUpperBoundaryForQx;

	EGS_RandomGenerator* m_pRandomNumberGenerator;

	unsigned int m_nNumberOfDeletedParticles;
};

#endif/*



*/
