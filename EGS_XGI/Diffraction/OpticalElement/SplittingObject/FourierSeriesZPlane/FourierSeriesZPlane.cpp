/*
###############################################################################
#
#   EGS_XGI FourierSeriesZPlane
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

#include "FourierSeriesZPlane.h"
#include "xgi_global_variables.h"
#include "egs_interface2.h"
#include "egs_base_geometry.h"
#include "egs_interpolator.h"


FourierSeriesZPlane::FourierSeriesZPlane(std::string i_sName, unsigned int i_nSplittingNumber, EGS_Float i_fPosition)
	:SplittingObject(i_sName, "FourierSeriesZPlane"),
	m_fPosition(i_fPosition),
	m_nNumberOfSecondaryPhotons(0),
	m_nNumberOfChargedParticles(0)
{

}

FourierSeriesZPlane::FourierSeriesZPlane(std::string i_sName)
	:SplittingObject(i_sName, "FourierSeriesZPlane"),
	m_nNumberOfSecondaryPhotons(0),
	m_nNumberOfChargedParticles(0)
{

}

FourierSeriesZPlane::FourierSeriesZPlane(std::string i_sName, EGS_Interpolator* i_i_gmfp, RefractiveIndexCalculator* i_pRefractiveIndexCalculator)
	:SplittingObject(i_sName, "FourierSeriesZPlane"),
	m_nNumberOfSecondaryPhotons(0),
	m_nNumberOfChargedParticles(0)
{
	i_gmfp = i_i_gmfp;
	m_bHavePrecalculatedCoefficients = false;
	m_pRefractiveIndexCalculator = i_pRefractiveIndexCalculator;
}

FourierSeriesZPlane::~FourierSeriesZPlane()
{

}

bool FourierSeriesZPlane::IsOnOpticalElement(EGS_Vector* i_pPosition)
{
	if(abs(i_pPosition->z - m_fPosition) < 1e-20)
	{
		return true;
	}
	else
	{
		return false;
	}
}




bool FourierSeriesZPlane::ApplyOpticalElementRule()
{
	m_nStackSizeBefore = the_stack->np - 1;
	EGS_Float fz_Before = the_stack->z[m_nStackSizeBefore];
	if (abs(fz_Before - m_fPosition) < 1e-10)
	{
		//On SplittingObject
		EGS_I32 nir_Before = the_stack->ir[m_nStackSizeBefore];
		int nCharge = the_stack->iq[m_nStackSizeBefore];
		EGS_Vector fDirection;
		fDirection.x = the_stack->u[m_nStackSizeBefore];
		fDirection.y = the_stack->v[m_nStackSizeBefore];
		fDirection.z = the_stack->w[m_nStackSizeBefore];
		fDirection.normalize();
		//only split:
		if(e_bPrimary[m_nStackSizeBefore] == true && nir_Before > 1 && nCharge == 0 && abs(fDirection.z) > 0.5 )
		{
			EGS_Float fEnergyBefore = the_stack->E[m_nStackSizeBefore];
			//Check if there are precalculated coefficients available
			if (m_bHavePrecalculatedCoefficients == false)
			{
				//calculate the coefficients
				StartRayTracing(fEnergyBefore);
			}
			//Get the relevant parameters from the_stack
			EGS_Float fPhase_Before = e_fPhase[m_nStackSizeBefore];
			EGS_Float fLogNorm_Before = e_fLogNorm[m_nStackSizeBefore];

			EGS_Float fwt_before = the_stack->wt[m_nStackSizeBefore];

			EGS_Float fx_Before = the_stack->x[m_nStackSizeBefore];
			EGS_Float fy_Before = the_stack->y[m_nStackSizeBefore];

			EGS_Float fdnear_Before = the_stack->dnear[m_nStackSizeBefore];

			EGS_I32 nlatch_Berfore = the_stack->latch[m_nStackSizeBefore];

			unsigned int nNumberOfCreatedPaths = 0;
			//Add all paths with relevant weight (Fourier coefficient) to the_stack
			if (m_bInclude_Coefficient_0 == true)
			{
				//Add the 0th Fourier coefficient path to the_stack
				e_fPhase[m_nStackSizeBefore] = fPhase_Before + m_fPhaseAmplitude_0;
				the_stack->wt[m_nStackSizeBefore] = fwt_before * m_fNormAmplitudeSquared_0 * m_fNormalizationConstant;

				nNumberOfCreatedPaths++;
			}
			if (m_bInclude_Coefficients_N == true)
			{
				//Add the all relevant Nth Fourier coefficient (N != 0) to the_stack (up to specified highest Fourier coefficient)
				for(unsigned int j = 0; j < m_fSin_over_npi_squared.size(); j++)
				{
					e_fPhase[m_nStackSizeBefore + nNumberOfCreatedPaths] = fPhase_Before - m_fQ_over_2pihbar[j] * fx_Before * m_fCurrentWaveLength + 0.5 * m_nIncludedFourierCoefficients[j] * m_fDutyCycle * m_fCurrentWaveLength + m_fPhaseAmplitude_N + m_fPhaseSin_over_npi[j] * m_fCurrentWaveLength;
					e_fLogNorm[m_nStackSizeBefore + nNumberOfCreatedPaths] = fLogNorm_Before;
					e_bPrimary[m_nStackSizeBefore + nNumberOfCreatedPaths] = true;

					the_stack->wt[m_nStackSizeBefore + nNumberOfCreatedPaths] = fwt_before * m_fSin_over_npi_squared[j] * m_fNormAmplitudeSquared_N * m_fNormalizationConstant;

					the_stack->E[m_nStackSizeBefore + nNumberOfCreatedPaths] = fEnergyBefore;

					EGS_Float fDirectionNewU = fDirection.x - m_fQ_over_2pihbar[j] * m_fCurrentWaveLength /*/ 2.0 / ec_fPi*/;
					the_stack->u[m_nStackSizeBefore + nNumberOfCreatedPaths] = fDirectionNewU;
					the_stack->v[m_nStackSizeBefore + nNumberOfCreatedPaths] = fDirection.y;
					the_stack->w[m_nStackSizeBefore + nNumberOfCreatedPaths] = sqrt(1.0 - fDirectionNewU * fDirectionNewU - fDirection.y * fDirection.y);

					the_stack->x[m_nStackSizeBefore + nNumberOfCreatedPaths] = fx_Before;
					the_stack->y[m_nStackSizeBefore + nNumberOfCreatedPaths] = fy_Before;
					the_stack->z[m_nStackSizeBefore + nNumberOfCreatedPaths] = fz_Before;

					the_stack->iq[m_nStackSizeBefore + nNumberOfCreatedPaths] = 0;
					the_stack->dnear[m_nStackSizeBefore + nNumberOfCreatedPaths] = fdnear_Before;

					the_stack->ir[m_nStackSizeBefore + nNumberOfCreatedPaths] = nir_Before;
					the_stack->latch[m_nStackSizeBefore + nNumberOfCreatedPaths] = nlatch_Berfore;

					nNumberOfCreatedPaths++;
				}
			}
			the_stack->np = m_nStackSizeBefore + nNumberOfCreatedPaths;

			//m_nSplittingNumber = nNumberOfCreatedPaths;
			m_nNumberOfTimesApplied++;
			m_nNumberOfGeneratedPaths += nNumberOfCreatedPaths;

			return true;
		}
		else if(nCharge == 0)
		{
			//only apply change of weight for photons
			//i.e. apply a uniform filter
			if (m_bIncludeAttenuationInGrating == true)
			{
				EGS_Float fDirectionCorrection = 1.0;
				if(m_bDirectionCorrectionGratingThickness == true && abs(the_stack->w[m_nStackSizeBefore]) > 0.0)
				{
					fDirectionCorrection = 1.0 / abs(the_stack->w[m_nStackSizeBefore]);
				}
				EGS_Float fNormTransmissionFunctionA = 1.0;
				EGS_Float fNormTransmissionFunctionB = 1.0;
				if (m_nMediumIndexSlabs >= 0)
				{
					fNormTransmissionFunctionA = exp(-1.0 * m_fHeightOfGratingSlabsInCM * fDirectionCorrection / i_gmfp[m_nMediumIndexSlabs].interpolate(log(the_stack->E[m_nStackSizeBefore])));
				}
				if (m_nMediumIndexGaps >= 0)
				{
					fNormTransmissionFunctionB = exp(-1.0 * m_fHeightOfGratingSlabsInCM * fDirectionCorrection / i_gmfp[m_nMediumIndexGaps].interpolate(log(the_stack->E[m_nStackSizeBefore])));
				}
				the_stack->wt[m_nStackSizeBefore] *= (m_fDutyCycle * fNormTransmissionFunctionA + (1.0 - m_fDutyCycle) * fNormTransmissionFunctionB);
				m_nNumberOfSecondaryPhotons++;
			}
			else
			{
				m_nNumberOfSecondaryPhotons++;
				the_stack->wt[m_nStackSizeBefore] *= (m_fDutyCycle * m_fNormTransmissionFunctionA + (1.0 - m_fDutyCycle) * m_fNormTransmissionFunctionB);
			}
			//m_nSplittingNumber = 1;
			return true;
		}
		else//electron or positron
		{
			//m_nSplittingNumber = 1;
			m_nNumberOfChargedParticles++;
			return true;
		}
	}
	else
	{
		return false;
	}
}


void FourierSeriesZPlane::ReportOpticalElement()
{
	std::cout << "--------------------------------" << std::endl;
	std::cout << "Splitting Object of Type: FourierSeriesZPlane" << std::endl;
	std::cout << "Name: " << m_sName << std::endl;
	std::cout << "Duty cycle: " << m_fDutyCycle << std::endl;
	std::cout << "Periodicity: " << m_fPeriodinCM << " cm" << std::endl;
	std::cout << "Highest order Fourier coefficient: " << m_nHighestFourierCoefficient << std::endl;
	std::cout << "Epsilon: " << m_fEpsilon << std::endl;
	std::cout << "Modelling of attenuation inside the grating is ";
	if (m_bIncludeAttenuationInGrating == true)
	{
		std::cout << "turned on." << std::endl;
		std::cout << "Height of the grating slabs: " << m_fHeightOfGratingSlabsInCM << " cm" << std::endl;
	}
	else
	{
		std::cout << "turned off." << std::endl;
	}
	std::cout << "Modelling directional dependence of grating thickness is ";
	if (m_bDirectionCorrectionGratingThickness == true)
	{
		std::cout << "turned on." << std::endl;
	}
	else
	{
		std::cout << "turned off." << std::endl;
	}
	std::cout << std::endl;
	std::cout << "--------------------------------" << std::endl;
}


int FourierSeriesZPlane::InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry)
{
	std::vector<EGS_Float> fPosition;
	std::vector<EGS_Float> fPhaseShiftAtDesignEnergy;//in units of pi
	std::vector<EGS_Float> fDutyCycle;
	std::vector<EGS_Float> fPeriodicityInCM;
	std::vector<int> nHighestFourierOrder;
	std::vector<EGS_Float> fDesignEnergy;//in MeV as in every other EGS fct

	//Get the corresponding values from EGS_Input
	int errPos = i_sInput->getInput("position", fPosition);
	int errPhaseShift = i_sInput->getInput("phase shift at design energy", fPhaseShiftAtDesignEnergy);
	int errDutyC = i_sInput->getInput("duty cycle", fDutyCycle);
	int errPeriod = i_sInput->getInput("periodicity", fPeriodicityInCM);
	int errFourierOrder = i_sInput->getInput("highest fourier order", nHighestFourierOrder);
	int errDesignEnergy = i_sInput->getInput("design energy", fDesignEnergy);


	if (errPos || errPhaseShift || errDutyC || errPeriod || errFourierOrder || errDesignEnergy)
	{
		std::cout << "FourierSeriesZPlane::InitOpticalElement: Error while reading parameters: 'position', 'phase shift at design energy', 'duty cycle', 'periodicity', 'highest fourier order', or 'design energy'." << std::endl;
		return 0;
	}
	else if (fPosition.size() == 1 && fPhaseShiftAtDesignEnergy.size() == 1 && fDutyCycle.size() == 1 && fPeriodicityInCM.size() == 1 && nHighestFourierOrder.size() == 1 && fDesignEnergy.size() == 1)
	{
		m_fPosition = fPosition[0];
		EGS_Float fDesignWaveLength = ec_fEnergyToWaveLength / fDesignEnergy[0];
		//Assuming quadratic dependence on wavelength for the real part of the refractive index
		m_PhaseConstantTransmissionFunctionA = ec_fPi * fPhaseShiftAtDesignEnergy[0] / fDesignWaveLength;
		m_fDutyCycle = fDutyCycle[0];
		m_fPeriodinCM = fPeriodicityInCM[0];
		m_nHighestFourierCoefficient = nHighestFourierOrder[0];
		if (m_nHighestFourierCoefficient < 0)
		{
			std::cout << "FourierSeriesZPlane::InitOpticalElement:Error: 'highest fourier order' must be larger than 0." << std::endl;
			return 0;
		}
	}
	else
	{
		std::cout << "FourierSeriesZPlane::InitOpticalElement: Unexpected format of parameters: 'position', 'phase shift at design energy', 'duty cycle', 'periodicity', 'highest fourier order', or 'design energy'." << std::endl;
		return 0;
	}

	std::vector<EGS_Float> fEpsilon;
	int errEpsilon = i_sInput->getInput("Epsilon", fEpsilon);
	if (errEpsilon)
	{
		std::cout << "FourierSeriesZPlane::InitOpticalElement: FourierSeriesZPlane::m_fEpsilon set to default value." << std::endl;
		m_fEpsilon = 0.0001;
	}
	else if (fEpsilon.size() == 1)
	{
		m_fEpsilon = fEpsilon[0];
	}
	else
	{
		std::cout << "FourierSeriesZPlane::InitOpticalElement: Unsupported format of 'epsilon' key" << std::endl;
		return 0;
	}

	//Define norm of the transmission function
	//2 options: setting as constant
	//            or enabling attenuation in the slabs
	//Set the flag for modelling attenuation
	int nBoolModelAttenuationInsideGrating;
	//If the include attenuation need the materials:
	std::string sMediumA;
	std::string sMediumB;
	EGS_Float fHeightOfGratingSlabs;
	//If don't include attenuation: need to define the norm of the transmission amplitude
	std::vector<EGS_Float> fNormOfTransmissionFunctionAtDesignEnergy;
	int nBoolDirectionCorrectionBasedOnInclommingDirection;

	int errAttenuation = i_sInput->getInput("model attenuation", nBoolModelAttenuationInsideGrating);
	int errMediumA = i_sInput->getInput("medium A", sMediumA);
	int errMediumB = i_sInput->getInput("medium B", sMediumB);
	int errHeight = i_sInput->getInput("grating thickness", fHeightOfGratingSlabs);
	int errNormTranmission = i_sInput->getInput("norm at design energy", fNormOfTransmissionFunctionAtDesignEnergy);
	int errDirectionCorrection = i_sInput->getInput("direction correction", nBoolDirectionCorrectionBasedOnInclommingDirection);


	if(errAttenuation)
	{
		std::cout << "FourierSeriesZPlane::InitOpticalElement: Assuming energy independent transmission function norm." << std::endl;
    m_bIncludeAttenuationInGrating = false;
	}
	else
	{
		if(nBoolModelAttenuationInsideGrating == 0)
		{
			m_bIncludeAttenuationInGrating = false;
		}
		else
		{
			m_bIncludeAttenuationInGrating = true;
		}
	}
	if (errDirectionCorrection)
	{
		std::cout << "FourierSeriesZPlane::InitOpticalElement: Error reading 'direction correction' (0->false, (!=0)->true)." << std::endl;
		return 0;
	}
	else
	{
		if (nBoolDirectionCorrectionBasedOnInclommingDirection == 0)
		{
			m_bDirectionCorrectionGratingThickness = false;
		}
		else //nBoolDirectionCorrectionBasedOnInclommingDirection == 1 or higher
		{
			m_bDirectionCorrectionGratingThickness = true;
		}
	}
	if(m_bIncludeAttenuationInGrating == false)
	{
		if (errNormTranmission)
    {
			std::cout << "FourierSeriesZPlane::InitOpticalElement: error reading 'norm at design energy' key." << std::endl;
      return 0;
    }
    else if (fNormOfTransmissionFunctionAtDesignEnergy.size() == 2)
    {
			m_fNormTransmissionFunctionA = fNormOfTransmissionFunctionAtDesignEnergy[0];
      m_fNormTransmissionFunctionB = fNormOfTransmissionFunctionAtDesignEnergy[1];
    }
    else
    {
			std::cout << "FourierSeriesZPlane::InitOpticalElement: Unexpected format of 'norm at design energy' key." << std::endl;
      return 0;
    }
	}
	else//m_bIncludeAttenuationInGrating == true
	{
		//Get paramaters for modelling energy dependent transmission function
		if(errMediumA || errMediumB || errHeight)
		{
			std::cout << "FourierSeriesZPlane::InitOpticalElement: error reading keys: 'medium A','medium B' or 'grating thickness'." << std::endl;
			return 0;
		}
		else
		{
			m_fHeightOfGratingSlabsInCM = fHeightOfGratingSlabs;
			m_nMediumIndexSlabs = i_pGeometry->getMediumIndex(sMediumA);
			std::cout << "FourierSeriesZPlane::InitOpticalElement:Slabs made of: " << sMediumA << "; " << m_nMediumIndexSlabs << std::endl;
			if (m_nMediumIndexSlabs == -1 && sMediumA.compare("VACUUM") != 0)
      {
				std::cout << "###############################################" << std::endl;
        std::cout << "FourierSeriesZPlane::InitOpticalElement:WARINING 'medium A' not defined in Geometry. Assume 'VACUUM' instead" << std::endl;
				std::cout << "###############################################" << std::endl;
      }
			m_nMediumIndexGaps = i_pGeometry->getMediumIndex(sMediumB);
			std::cout << "FourierSeriesZPlane::InitOpticalElement:Gaps filled with " << sMediumB << "; " << m_nMediumIndexGaps << std::endl;
			if (m_nMediumIndexGaps == -1 && sMediumB.compare("VACUUM") != 0)
      {
				std::cout << "###############################################" << std::endl;
      	std::cout << "FourierSeriesZPlane::InitOpticalElement:WARNING 'medium B' not defined in Geometry. Assume 'VACUUM' instead" << std::endl;
				std::cout << "###############################################" << std::endl;
      }
			if (m_nMediumIndexSlabs == -1 && m_nMediumIndexGaps == -1)
      {
				std::cout << "###############################################" << std::endl;
      	std::cout << "FourierSeriesZPlane::InitOpticalElement:WARNING: Grating attenuation modelling is turned off because both sections are set to 'VACUUM'." << std::endl;
				std::cout << "###############################################" << std::endl;
        m_bIncludeAttenuationInGrating = false;
      }
      if (m_nMediumIndexSlabs == -1)
      {
        std::cout <<"FourierSeriesZPlane::InitOpticalElement: set norm of transmission function 'A' inside slabs const 1.0." << std::endl;
        m_fNormTransmissionFunctionA = 1.0;
      }
      if (m_nMediumIndexGaps == -1)
      {
        std::cout << "FourierSeriesZPlane::InitOpticalElement: set norm of transmission function 'B' between slabs const 1.0." << std::endl;
        m_fNormTransmissionFunctionB = 1.0;
      }
		}
	}

	//Check which Fourier coefficients have to be included in the splitting
	//Init sin(n pi a/(a+b))/(n pi)
	m_fNormalizationConstant_N = 0.0;
	for (int nFourierCoefficientCounter = 1; nFourierCoefficientCounter <= m_nHighestFourierCoefficient; nFourierCoefficientCounter++)
	{
		EGS_Float fSin_N_Pi_DC = sin(nFourierCoefficientCounter * ec_fPi * m_fDutyCycle);
		//If the sin^2 > m_fEpsilon take into account the corresponding Fourier coefficient in splitting.
		if (fSin_N_Pi_DC * fSin_N_Pi_DC > m_fEpsilon)
		{
			EGS_Float fEnergyIndependentPartOfAmplitudef = fSin_N_Pi_DC / ec_fPi / nFourierCoefficientCounter;
			EGS_Float fQ_over_2pihbar = /*2 * ec_fPi **/ nFourierCoefficientCounter / m_fPeriodinCM;
			m_fSin_over_npi_squared.push_back(fEnergyIndependentPartOfAmplitudef * fEnergyIndependentPartOfAmplitudef);
			m_fQ_over_2pihbar.push_back(fQ_over_2pihbar);
			m_nIncludedFourierCoefficients.push_back(nFourierCoefficientCounter);

			m_fSin_over_npi_squared.push_back(fEnergyIndependentPartOfAmplitudef * fEnergyIndependentPartOfAmplitudef);
			m_fQ_over_2pihbar.push_back((-1.0) * fQ_over_2pihbar);
			m_nIncludedFourierCoefficients.push_back(-nFourierCoefficientCounter);

			m_fNormalizationConstant_N += 2.0 * fEnergyIndependentPartOfAmplitudef * fEnergyIndependentPartOfAmplitudef;

			if (fEnergyIndependentPartOfAmplitudef < 0)
      {
      	m_fPhaseSin_over_npi.push_back(0.5);//fPi/2/pi ;
				m_fPhaseSin_over_npi.push_back(0.5);//Phase of +-N are the same
      }
      else
      {
      	m_fPhaseSin_over_npi.push_back(0.0);
				m_fPhaseSin_over_npi.push_back(0.0);
      }
		}
	}
	return 1;
}


void FourierSeriesZPlane::EndRayTracing()
{
	m_bHavePrecalculatedCoefficients = false;
}

//Calculate energy dependent part of Fourier coefficients
void FourierSeriesZPlane::StartRayTracing(EGS_Float& i_fEnergy)
{
	//If no thickness correction is made "lock" the calculation of the energy dependent part of the coefficients
		//otherwise the correction has to be applied for each path.

	//The correction factor for the grating thickness
	EGS_Float fDirectionCorrection = 1.0;
	//Constant phase for each splitting
	EGS_Float fPhaseOfPrefactor = 0.0;
	if (m_bDirectionCorrectionGratingThickness == false)
	{
		m_bHavePrecalculatedCoefficients = true;
	}
	else
	{
		int nStackAddress = the_stack->np - 1;
		if (the_stack->w[nStackAddress] > 0.0)
		{
			fDirectionCorrection = 1.0 / the_stack->w[nStackAddress];
			fPhaseOfPrefactor = (1.0 - m_pRefractiveIndexCalculator->GetRefractiveIndexDelta(m_nMediumIndexGaps, the_stack->E[nStackAddress])) * m_fHeightOfGratingSlabsInCM * fDirectionCorrection;
		}
		else
		{
			std::cout << "FourierSeriesZPlane::StartRayTracing: WARNING: Path parallel to grating." << std::endl;
		}
	}

	//The wave length
	m_fCurrentWaveLength = ec_fEnergyToWaveLength / the_stack->E[m_nStackSizeBefore];

	EGS_Float fPhiA = m_PhaseConstantTransmissionFunctionA * m_fCurrentWaveLength * fDirectionCorrection;

	if (m_bIncludeAttenuationInGrating == true)
	{
		if (m_nMediumIndexSlabs >= 0)
		{
			m_fNormTransmissionFunctionA = exp(-1.0 * m_fHeightOfGratingSlabsInCM * fDirectionCorrection / i_gmfp[m_nMediumIndexSlabs].interpolate(log(i_fEnergy)) / 2.0);
		}
		if (m_nMediumIndexGaps >= 0)
		{
			m_fNormTransmissionFunctionB = exp(-1.0 * m_fHeightOfGratingSlabsInCM * fDirectionCorrection / i_gmfp[m_nMediumIndexGaps].interpolate(log(i_fEnergy)) / 2.0);
		}
	}

	//Calculate energy dependent part of the Fourier coefficients: |\tau_a - \tau_b|^2 and |a\tau_a + b\tau_b|^2/|a+b|^2
	//where \tau_a = m_fNormTransmissionFunctionA * exp(I fPhiA)
	//and \tau_b = m_fNormTransmissionFunctionA * exp(I * 0)
	EGS_Float fEnergyDependentPartOfFC_N_Real = m_fNormTransmissionFunctionA * cos(fPhiA) - m_fNormTransmissionFunctionB;
	EGS_Float fEnergyDependentPartOfFC_N_Imag = m_fNormTransmissionFunctionA * sin(fPhiA);
	EGS_Float fEnergyDependentPartOfFC_0_Real = m_fNormTransmissionFunctionA * cos(fPhiA) * m_fDutyCycle + m_fNormTransmissionFunctionB * (1.0 - m_fDutyCycle);
	EGS_Float fEnergyDependentPartOfFC_0_Imag = m_fNormTransmissionFunctionA * sin(fPhiA) * m_fDutyCycle;

	m_fNormAmplitudeSquared_N = fEnergyDependentPartOfFC_N_Real * fEnergyDependentPartOfFC_N_Real + fEnergyDependentPartOfFC_N_Imag * fEnergyDependentPartOfFC_N_Imag;
	m_fNormAmplitudeSquared_0 = fEnergyDependentPartOfFC_0_Real * fEnergyDependentPartOfFC_0_Real + fEnergyDependentPartOfFC_0_Imag * fEnergyDependentPartOfFC_0_Imag;

	EGS_Float fTargetNormalization = m_fDutyCycle * m_fNormTransmissionFunctionA * m_fNormTransmissionFunctionA + (1.0 - m_fDutyCycle) * m_fNormTransmissionFunctionB * m_fNormTransmissionFunctionB;

	if (m_fNormAmplitudeSquared_N > 0 || m_fNormAmplitudeSquared_0 > 0)
	{
		if (m_fNormAmplitudeSquared_N > m_fNormAmplitudeSquared_0)
		{
			if (m_fNormAmplitudeSquared_0 / m_fNormAmplitudeSquared_N < m_fEpsilon)
			{
				//Ignore coefficient 0
				m_bInclude_Coefficient_0 = false;
				m_bInclude_Coefficients_N = true;
				m_fNormalizationConstant = fTargetNormalization / (m_fNormalizationConstant_N * m_fNormAmplitudeSquared_N);
			}
			else
			{
				//include all coefficients
				m_bInclude_Coefficient_0 = true;
				m_bInclude_Coefficients_N = true;
				m_fNormalizationConstant = fTargetNormalization / (m_fNormAmplitudeSquared_0 + m_fNormalizationConstant_N * m_fNormAmplitudeSquared_N);
			}
		}
		else
		{
			if (m_fNormAmplitudeSquared_N / m_fNormAmplitudeSquared_0 < m_fEpsilon)
			{
				m_bInclude_Coefficient_0 = true;
				m_bInclude_Coefficients_N = false;
				m_fNormalizationConstant = fTargetNormalization / (m_fNormAmplitudeSquared_0);
			}
			else
			{
				m_bInclude_Coefficient_0 = true;
				m_bInclude_Coefficients_N = true;
				m_fNormalizationConstant = fTargetNormalization / (m_fNormAmplitudeSquared_0 + m_fNormalizationConstant_N * m_fNormAmplitudeSquared_N);
			}
		}
		if (m_bInclude_Coefficient_0 == true)
		{
			m_fPhaseAmplitude_0 = GetPhaseOfComplexNumber(fEnergyDependentPartOfFC_0_Real, fEnergyDependentPartOfFC_0_Imag) * m_fCurrentWaveLength / 2.0 / ec_fPi + fPhaseOfPrefactor;
		}
		else
		{
			m_fPhaseAmplitude_0 = 0.0;
		}
		if (m_bInclude_Coefficients_N == true)
		{
			m_fPhaseAmplitude_N = GetPhaseOfComplexNumber(fEnergyDependentPartOfFC_N_Real, fEnergyDependentPartOfFC_N_Imag) * m_fCurrentWaveLength / 2.0 / ec_fPi + fPhaseOfPrefactor;
		}
		else
		{
			m_fPhaseAmplitude_N = 0.0;
		}
	}
	else
	{
		std::cout << "FourierSeriesZPlane::StartRayTracing: Error unable to calculate Fourier coeficients" << std::endl;
		m_bInclude_Coefficient_0 = true;
		m_fPhaseAmplitude_0 = 0.0;
		m_bInclude_Coefficients_N = false;
		m_fPhaseAmplitude_N = 0.0;
	}
}

void FourierSeriesZPlane::PrintSummary()
{
	std::cout << "--------------------------------" << std::endl;
	std::cout << "Splitting Object of Type: " << m_sType << std::endl;
	std::cout << "Name: " << m_sName << std::endl;
	std::cout << "Number of rule applied: " << m_nNumberOfTimesApplied << std::endl;
	std::cout << "Number of paths created: " << m_nNumberOfGeneratedPaths << std::endl;
	std::cout << "Number of secondary photons: "  << m_nNumberOfSecondaryPhotons << std::endl;
	std::cout << "Number of charged particles: "  << m_nNumberOfChargedParticles << std::endl;
	std::cout << "--------------------------------" << std::endl;
}


bool FourierSeriesZPlane::PermissionToResetPointer(int i_nStackSize)
{
	if(m_nStackSizeBefore == i_nStackSize)
	{
		return true;
	}
	else if(m_nStackSizeBefore >= i_nStackSize)
	{
		std::cout << "SplittingObject::PermissionToResetPointer: "<< m_sName << ": warning stack size too low: " << m_nStackSizeBefore << " > " << i_nStackSize << std::endl;
		//std::cout << "m_nSplittingNumber = " << m_nSplittingNumber << std::endl;
		std::cout << "latch: " << the_stack->latch[i_nStackSize] << std::endl;
		std::cout << "x = " << the_stack->x[i_nStackSize] << std::endl;
		std::cout << "y = " << the_stack->y[i_nStackSize] << std::endl;
		std::cout << "z = " << the_stack->z[i_nStackSize] << std::endl;
		std::cout << "u = " << the_stack->u[i_nStackSize] << std::endl;
		std::cout << "v = " << the_stack->v[i_nStackSize] << std::endl;
		std::cout << "w = " << the_stack->w[i_nStackSize] << std::endl;
		std::cout << "wt = " << the_stack->wt[i_nStackSize] << std::endl;
		std::cout << "E = " << the_stack->E[i_nStackSize] << std::endl;
		std::cout << "iq = " << the_stack->iq[i_nStackSize] << std::endl;
		std::cout << "ir = " << the_stack->ir[i_nStackSize] - 2 << std::endl;
		return true;
	}
	else
	{
		return false;
	}
}
