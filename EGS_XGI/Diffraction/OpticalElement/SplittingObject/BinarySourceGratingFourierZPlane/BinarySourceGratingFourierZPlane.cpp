/*
###############################################################################
#
#   EGS_XGI BinarySourceGratingFourierZPlane
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
#include "BinarySourceGratingFourierZPlane.h"

#include "xgi_global_variables.h"

#include "egs_interface2.h"

BinarySourceGratingFourierZPlane::BinarySourceGratingFourierZPlane(std::string i_sName, EGS_RandomGenerator* i_pRandomNumberGenerator)
	:SplittingObject(i_sName, "BinarySourceGratingFourierZPlane"),
	m_pRandomNumberGenerator(i_pRandomNumberGenerator)
{
	m_nNumberOfDeletedParticles = 0;
}

BinarySourceGratingFourierZPlane::BinarySourceGratingFourierZPlane(std::string i_sName)
	:SplittingObject(i_sName, "BinarySourceGratingFourierZPlane")
{
	m_nNumberOfDeletedParticles = 0;
}

BinarySourceGratingFourierZPlane::~BinarySourceGratingFourierZPlane()
{

}

int BinarySourceGratingFourierZPlane::InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry)
{
	std::vector<EGS_Float> fPosition;
	std::vector<EGS_Float> fTransmissionFunctionNorm;//RE + i*IM
	std::vector<EGS_Float> fHoleSizeInCM;
	std::vector<EGS_Float> fPeriodicityInCM;
	std::vector<EGS_Float> fHighestZeroToConsider;
	std::vector<int> nNumberOfPathsToCreate;

	int errPos = i_sInput->getInput("position", fPosition);
	int errTransA = i_sInput->getInput("transmission function norm a", fTransmissionFunctionNorm);
	int errHSize = i_sInput->getInput("slit width", fHoleSizeInCM);
	int errPeriod = i_sInput->getInput("periodicity", fPeriodicityInCM);
	int errZero = i_sInput->getInput("highest zero", fHighestZeroToConsider);//up to which zero of the sinc the Qx are considered (can choose fractions...)
	int errClones = i_sInput->getInput("number of paths", nNumberOfPathsToCreate);

	if(errPos || errTransA || errHSize || errPeriod || errZero || errClones)
	{
		std::cout << "BinarySourceGratingFourierZPlane::InitOpticalElement: Error while reading parameters" << std::endl;
		return 0;
	}
	else if(fPosition.size() == 1 && fHoleSizeInCM.size() == 1 && fPeriodicityInCM.size() == 1 && fHighestZeroToConsider.size() == 1 && nNumberOfPathsToCreate.size() == 1)
	{
		m_fPosition = fPosition[0];

		m_fTransmissionNorm = fTransmissionFunctionNorm[0];
		m_fHoleSizeInCM = fHoleSizeInCM[0];
		m_fPeriodicityInCM = fPeriodicityInCM[0];

		m_fUpperBoundaryForQx = 2.0 * ec_fPi * fHighestZeroToConsider[0] / m_fHoleSizeInCM;
		m_nSplittingNumber = 2 * nNumberOfPathsToCreate[0];

		return 1;
	}
	else
	{
		std::cout << "BinarySourceGratingFourierZPlane::InitOpticalElement: Error while reading parameters" << std::endl;
		return 0;
	}
}

bool BinarySourceGratingFourierZPlane::IsOnOpticalElement (EGS_Vector* i_pPosition)
{
	if(abs(i_pPosition->z - m_fPosition) < 1e-10)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool BinarySourceGratingFourierZPlane::ApplyOpticalElementRule()
{
	m_nStackSizeBefore = the_stack->np - 1;
	EGS_Float fZPosition = the_stack->z[m_nStackSizeBefore];
	if(abs(fZPosition - m_fPosition) > 1e-10)
	{
		return false;
	}
	else
	{
		EGS_Float fXPosition = the_stack->x[m_nStackSizeBefore];
		int nGratingSection = floor(fXPosition / m_fPeriodicityInCM);
		if(fXPosition - double(nGratingSection) * m_fPeriodicityInCM > m_fHoleSizeInCM)
		{
			//The particle hits a perfectly absorbing grating slab
			//Make sure the particle gets deleted
			the_stack->z[m_nStackSizeBefore] = 1e30;
			the_stack->ir[m_nStackSizeBefore] = -1;
			m_nNumberOfDeletedParticles++;
			return false;
		}
		else
		{
			//pass through a slit
			EGS_I32 nRegionIndexBefore = the_stack->ir[m_nStackSizeBefore];
			EGS_Vector fDirection;
			fDirection.x = the_stack->u[m_nStackSizeBefore];
			fDirection.y = the_stack->v[m_nStackSizeBefore];
			fDirection.z = the_stack->w[m_nStackSizeBefore];
			fDirection.normalize();
			if(e_bPrimary[m_nStackSizeBefore] == true && nRegionIndexBefore > 1 && the_stack->iq[m_nStackSizeBefore] == 0 && abs(fDirection.z) > 0.5)
			{
				EGS_Float fWeightBeforePathSplitting = the_stack->wt[m_nStackSizeBefore];
				EGS_Float fPhaseBeforePathSplitting = e_fPhase[m_nStackSizeBefore];
				EGS_Float fNormBeforeSpitting = e_fLogNorm[m_nStackSizeBefore];

				//EGS_Float u_Before = the_stack-> u[m_nStackSizeBefore];
				//EGS_Float v_Before = the_stack-> v[m_nStackSizeBefore];

				EGS_Float fEnergyBefore = the_stack->E[m_nStackSizeBefore];
				EGS_Float fWavelength_over_2pi = ec_fEnergyToWaveLength / (fEnergyBefore * 2.0 * ec_fPi);

				EGS_Float fYPosition = the_stack->y[m_nStackSizeBefore];

				EGS_Float dnearBefore = the_stack->dnear[m_nStackSizeBefore];
				//EGS_I32 nRegionIndexBefore = the_stack->ir[m_nStackSizeBefore];
				EGS_I32 nLatchBerfore = the_stack->latch[m_nStackSizeBefore];

				EGS_Float fNormalizationConstantProbabilities = 0.0;

				for(int nPathCounter = 0; nPathCounter < m_nSplittingNumber / 2; nPathCounter++)
				{
					EGS_Float fDirectionChangeX = m_pRandomNumberGenerator->getUniform() * m_fUpperBoundaryForQx;
					EGS_Float fSinc = 0.0;
					if(abs(fDirectionChangeX) < 1e-10)
					{
						//Use 2nd order Taylor expansion to compute $\tilde\tau(Q_x)$
						fSinc = m_fHoleSizeInCM / 2.0 - m_fHoleSizeInCM * m_fHoleSizeInCM * m_fHoleSizeInCM * fDirectionChangeX * fDirectionChangeX / 48.0;
					}
					else
					{
						fSinc = sin(fDirectionChangeX * m_fHoleSizeInCM / 2.0) / fDirectionChangeX;
					}
					EGS_Float fProbability =  fSinc * fSinc * m_fTransmissionNorm;
					fNormalizationConstantProbabilities += 2.0 * fSinc * fSinc;
					EGS_Float fPhase = fWavelength_over_2pi * fDirectionChangeX * (EGS_Float(nGratingSection) * m_fPeriodicityInCM  + m_fHoleSizeInCM / 2.0 - fXPosition);

					e_fPhase[m_nStackSizeBefore + nPathCounter] = fPhaseBeforePathSplitting + fPhase;
					e_fLogNorm[m_nStackSizeBefore + nPathCounter] = fNormBeforeSpitting;
					e_bPrimary[m_nStackSizeBefore + nPathCounter] = true;

					the_stack->E[m_nStackSizeBefore + nPathCounter] 	= fEnergyBefore;

					the_stack->x[m_nStackSizeBefore + nPathCounter] 	= fXPosition;
					the_stack->y[m_nStackSizeBefore + nPathCounter] 	= fYPosition;
					the_stack->z[m_nStackSizeBefore + nPathCounter] 	= fZPosition;

					the_stack->u[m_nStackSizeBefore + nPathCounter] 	= fDirection.x - fDirectionChangeX * fWavelength_over_2pi;
					the_stack->v[m_nStackSizeBefore + nPathCounter] 	= fDirection.y;
					the_stack->w[m_nStackSizeBefore + nPathCounter] 	= sqrt(1.0 - the_stack->u[m_nStackSizeBefore + nPathCounter] * the_stack->u[m_nStackSizeBefore + nPathCounter] - fDirection.y * fDirection.y);

					the_stack->dnear[m_nStackSizeBefore + nPathCounter] = dnearBefore;

					the_stack->wt[m_nStackSizeBefore + nPathCounter] 	= fWeightBeforePathSplitting * fProbability;

					the_stack->iq[m_nStackSizeBefore + nPathCounter] 	= (EGS_I32)0;
					the_stack->ir[m_nStackSizeBefore + nPathCounter] 	= nRegionIndexBefore;
					the_stack->latch[m_nStackSizeBefore + nPathCounter] = nLatchBerfore;


					//minus Q
					e_fPhase[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] = fPhaseBeforePathSplitting - fPhase;
					e_fLogNorm[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] = fNormBeforeSpitting;
					e_bPrimary[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] = true;

					the_stack->E[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] 	= fEnergyBefore;

					the_stack->x[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] 	= fXPosition;
					the_stack->y[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] 	= fYPosition;
					the_stack->z[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] 	= fZPosition;

					the_stack->u[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] 	= fDirection.x + fDirectionChangeX * fWavelength_over_2pi;
					the_stack->v[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] 	= fDirection.y;
					the_stack->w[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] 	= sqrt(1.0 - the_stack->u[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] * the_stack->u[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] - fDirection.y * fDirection.y);

					the_stack->dnear[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] = dnearBefore;

					the_stack->wt[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] = fWeightBeforePathSplitting * fProbability;

					the_stack->iq[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2]	= (EGS_I32)0;
					the_stack->ir[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2]	= nRegionIndexBefore;
					the_stack->latch[m_nStackSizeBefore + nPathCounter + m_nSplittingNumber / 2] = nLatchBerfore;
				}

				if(fNormalizationConstantProbabilities <= 0.0)
				{
					std::cout << "Warning fNormalizationConstantProbabilities = "  << fNormalizationConstantProbabilities << std::endl;
				}
				else
				{
					fNormalizationConstantProbabilities = 1.0 / fNormalizationConstantProbabilities;
					for(int nPathCounter = 0; nPathCounter < m_nSplittingNumber; nPathCounter++)
					{
						the_stack->wt[m_nStackSizeBefore + nPathCounter] *= fNormalizationConstantProbabilities;
					}
				}
				the_stack->np += m_nSplittingNumber -1;
				//the_stack->npold = the_stack->np;
				m_nNumberOfGeneratedPaths += m_nSplittingNumber -1;

				m_nNumberOfTimesApplied++;

				//EGS_Float fCheckNormalization = 0.0;
				return true;
			}
			else
			{
				//A secondary particle went through a slit
				//or a stack entry with ireg = 1 (outside) is processed
				return true;
			}
		}
	}
}


void BinarySourceGratingFourierZPlane::ReportOpticalElement()
{
	std::cout << "--------------------------------" << std::endl;
	std::cout << "Splitting Object of Type: BinarySourceGratingFourierZPlane" << std::endl;
	std::cout << "Name: " << m_sName << std::endl;
	std::cout << "Splitting number = " << m_nSplittingNumber << std::endl;
	std::cout << "Periodicity = " << m_fPeriodicityInCM << std::endl;
	std::cout << "Slit width = " << m_fHoleSizeInCM << std::endl;
	std::cout << "Transmission funciton norm = " << m_fTransmissionNorm << std::endl;
	std::cout << "--------------------------------" << std::endl;
}


void BinarySourceGratingFourierZPlane::PrintSummary()
{
	std::cout << "--------------------------------" << std::endl;
	std::cout << "Splitting Object of Type: " << m_sType << std::endl;
	std::cout << "Name: " << m_sName << std::endl;
	std::cout << "Number of splittings: " << m_nNumberOfTimesApplied << std::endl;
	std::cout << "Number of paths created: " << m_nNumberOfGeneratedPaths << std::endl;
	std::cout << "Number of absorbed stack entries: " << m_nNumberOfDeletedParticles << std::endl;
	std::cout << "--------------------------------" << std::endl;
}
