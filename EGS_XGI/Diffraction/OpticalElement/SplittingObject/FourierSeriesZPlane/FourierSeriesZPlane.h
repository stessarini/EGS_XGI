/*
###############################################################################
#
#   EGS_XGI FourierSeriesZPlane header
#   Model phase gratings (orthogonal to z-axis) with Fourier series
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

#ifndef _FOURIERSERIESZPLANE_H_
#define _FOURIERSERIESZPLANE_H_

#include <fstream>
#include <vector>


#include "SplittingObject.h"
#include "RefractiveIndexCalculator.h"

/*A class to simulate phase gratings.
  Splits incoming paths into paths according to the Fourier series of the
  transmission function. The argument of the Fourier coefficients gives the
  change of direction of the path.

  Has modes selected by the boolean member variables:
  - 'm_bIncludeAttenuationInGrating'
	->Switches on $Exp[-\mu * d / 2]$-terms in the transmission function,
	   where $\mu$ denotes the attenuation coefficient and d the grating
	   thickness.
  - 'm_bDirectionCorrectionGratingThickness'
	->Switches on a correction factor 'fDirectionCorrection' for the grating
	  thickness. 'fDirectionCorrection' = $1.0 / (\hat{p}\hat{n})$, where $\hat{p}$
	  denotes the direction of the momentum vector and $\hat{n}$ the normal
	  vector on the grating, both with unit length :
	  $\left|\hat{p}\right| = \left|\hat{n}\right| = 1$. In the case of a
	  flat grating 'fDirectionCorrection' = 1.0/the_stack->w[.].
  The two functionalities are independent of each other.
  Energy dependence of the phase shift is aways taken into account.

  The Fourier series is split into two parts an energy independent part
  $sin(n\pi 'm_fDutyCycle')/(n\pi)$ done in InitOpticalElement, an energy dependent
  part $\tau_a -\tau_b$ (done in StartRayTracing()) and a spatially dependent phase
  shift $Exp[-iQx_1]$, where $x_1$ is the x-position in the grating plane, calculated
  in ApplyOpticalElementRule().
  If 'm_bDirectionCorrectionGratingThickness' is set to false the Fourier series only
  has to be calculated once for each "incoherent history" (i.e. indistinguishable
  alternatives don't require a recalculation of the Fourier series).
  */
class FourierSeriesZPlane : public SplittingObject
{
public:
	FourierSeriesZPlane(std::string i_sName, unsigned int i_nSplittingNumber, EGS_Float i_fPosition);
	FourierSeriesZPlane(std::string i_sName);
	FourierSeriesZPlane(std::string i_sName, EGS_Interpolator* i_i_gmfp, RefractiveIndexCalculator* i_pRefractiveIndexCalculator);
	~FourierSeriesZPlane();

	virtual bool IsOnOpticalElement (EGS_Vector* i_pPosition);
	virtual bool ApplyOpticalElementRule();
	virtual void ReportOpticalElement();
	virtual void PrintSummary();

	/*Some part of the Fourier coefficients needed for the splitting is energy dependent.
		  This part can be reused in all the splittings concerning the same history. This function
		  informs the grating that next time a splitting occurs the coefficients have to be
		  recalculated.*/
	virtual void EndRayTracing();

	/* When first asked to split the incoming path: Calculate the energy dependent part of the Fourier coefficients
		   to enable faster transport. Once done set a flag that tells the grating to reuse
		   the calculated coefficients*/
	virtual void StartRayTracing(EGS_Float& i_fEnergy);

	/*Setup and calculate the history independent part of the Fourier coefficients
	  i.e. the sin(n \pi a/(a+b))/(n \pi)*/
	virtual int InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry);

	bool PermissionToResetPointer(int i_nStackSize);
private:
	EGS_Float m_fPosition;

	//Grating parameters
	EGS_Float m_fDutyCycle;
	EGS_Float m_fPeriodinCM;
	EGS_Float m_fPhaseShiftAtDesignEnergyInUnitsOfPi;
	//Assume quadratic dependence of phase shift
	//Hence use the following constants to determine the phase shift for different energies
	/*It is assumed that the input specifies the relative phase shift of the grating sections.
	  (additional phase shift in section a)*/
	EGS_Float m_PhaseConstantTransmissionFunctionA;

	//In case attenuation in the slabs is taken into account
	//Assume that the gaps are filled with air.
	//Assume that air attenuation can be neglected.
	bool m_bIncludeAttenuationInGrating;
	EGS_Float m_fHeightOfGratingSlabsInCM; //How high in z-direction
	int m_nMediumIndexSlabs; //Medium inside slab sections 'a'
	int m_nMediumIndexGaps; //Media in which the grating is embedded 'b'
	EGS_Interpolator* i_gmfp;//EGS interpolator for attenuation coefficients
	/*In case one is not taking into account attenuation inside grating slabs,
	set a default value for the transmission function.
	In that case the Norm of the transmission function is saved in the following
	two members*/
	EGS_Float m_fNormTransmissionFunctionA;
	EGS_Float m_fNormTransmissionFunctionB;
	EGS_Float m_fCurrentWaveLength;

	//Splitting related members
	//Constant for the whole transport
	std::vector<EGS_Float> m_fSin_over_npi_squared;
	std::vector<EGS_Float> m_fPhaseSin_over_npi;
	std::vector<EGS_Float> m_fQ_over_2pihbar;
	std::vector<int> m_nIncludedFourierCoefficients;
	EGS_Float m_fNormalizationConstant_N;
	EGS_Float m_fNormalizationConstant;

	//Constant for each MC history
	EGS_Float m_fNormAmplitudeSquared_N;
	EGS_Float m_fPhaseAmplitude_N;
	EGS_Float m_fNormAmplitudeSquared_0;
	EGS_Float m_fPhaseAmplitude_0;
	//Flags for evaluation
	bool m_bInclude_Coefficients_N;
	bool m_bInclude_Coefficient_0;
	bool m_bHavePrecalculatedCoefficients;

	//Flag for dependence of transmission function on direction of incoming path
	bool m_bDirectionCorrectionGratingThickness;

	//Numerical constant
	EGS_Float m_fEpsilon;

	unsigned int m_nNumberOfSecondaryPhotons;
	unsigned int m_nNumberOfChargedParticles;



	int m_nHighestFourierCoefficient;
  RefractiveIndexCalculator* m_pRefractiveIndexCalculator;
};

#endif
