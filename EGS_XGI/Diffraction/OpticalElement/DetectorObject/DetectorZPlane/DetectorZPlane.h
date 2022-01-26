/*
###############################################################################
#
#   EGS_XGI DetectorZPlane header
#   A basic detector class (orthogonal to z-axis).
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

#ifndef _DETECTORZPLANE_H_
#define _DETECTORZPLANE_H_

#include "OpticalElement.h"
#include <vector>
#include <string>
#include <fstream>


class DetectorZPlane : public OpticalElement
{
public:
	DetectorZPlane();

	DetectorZPlane(unsigned int i_nNxTotalSignal, unsigned int i_nNyTotalSignal, EGS_Vector& i_fPosition, EGS_Float i_fWidth, EGS_Float i_fHeight, bool i_bWirteToBinary, bool i_bWriteToPGM, std::string i_sFileName, std::string i_sDetectorName);

	DetectorZPlane(std::string i_sName, std::string i_sInputFileName);

	~DetectorZPlane();

	virtual bool IsOnOpticalElement (EGS_Vector* i_pPosition);
	virtual bool ApplyOpticalElementRule();
	virtual void PrintSummary();
	virtual void ReportOpticalElement();

	virtual void EndRayTracing();

	void WriteOnFile(long int i_nLastCase);
	void CreateImage();
	void WriteToFileNumberOfPaths();
	virtual void ScoreSecondaryParticle();
	virtual bool PermissionToResetPointer(int i_nStackSize);
	virtual int InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry);

	virtual void EndSimulation(int i_nNumberOfHistories);
protected:

	void ResetTotalSignal();

	std::vector<double>* m_pTotalSignal;
	std::vector<double>* m_pHistorySignalRealPart;
	std::vector<double>* m_pHistorySignalImaginaryPart;

	unsigned int m_nNxTotalSignal;
	unsigned int m_nNyTotalSignal;

		EGS_Vector m_fPosition;
	EGS_Float m_fWidth;
	EGS_Float m_fHeight;

	EGS_Float m_f_TotalSignal_DeltaX; //Bin width
	EGS_Float m_f_Totalsignal_DeltaY;


	bool m_bWirteToBinary;
	bool m_bWriteToPGM;

	std::string m_sFileName;

	bool m_bCountNumberOfPaths;
	std::vector<unsigned int>* m_pNumberOfPaths;

	std::string m_sPath;
	bool m_bAllowedToOverWrite;

	unsigned int m_nNumberOfSecondaryPhotons;

	bool m_bGenerateSeparateBinaryForSecondarySignal;
	bool m_bDoHistoryWiseNormalization;
	EGS_Float m_fTotalIncommingStatWeightOfPrimaryPaths;
};

#endif
