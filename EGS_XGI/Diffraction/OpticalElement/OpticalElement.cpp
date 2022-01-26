#include "OpticalElement.h"
#include "xgi_global_variables.h"
#include "egs_interpolator.h"
#include "egs_base_geometry.h"
OpticalElement::OpticalElement()
	:m_nNumberOfTimesApplied(0)
{
	m_nStackSizeBefore = -1;
}

OpticalElement::~OpticalElement()
{

}

bool OpticalElement::IsOnOpticalElement (EGS_Vector* i_pPosition)
{
	return false;
}

bool OpticalElement::ApplyOpticalElementRule()
{
	return true;
}

void OpticalElement::PrintSummary()
{

}

void OpticalElement::ReportOpticalElement()
{

}

EGS_Float OpticalElement::GetNumberOfTimesAppliedOpticalElemetRule()
{
	return m_nNumberOfTimesApplied;
}

void OpticalElement::EndRayTracing()
{

}

void OpticalElement::StartRayTracing(EGS_Float& i_fEnergy)
{

}

std::string OpticalElement::GetName()
{
	return m_sName;
}

void OpticalElement::ScoreSecondaryParticle()
{

}

int OpticalElement::InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry)
{
	return 0;
}

double OpticalElement::GetPhaseOfComplexNumber(EGS_Float i_fRealPart, EGS_Float i_fImaginaryPart)
{
	if(abs(i_fRealPart) <= 1e-20)
	{
		if(abs(i_fImaginaryPart) <= 1e-20)
		{
			return 0.0;
		}
		else if(i_fImaginaryPart >= 0.0)
		{
			return ec_fPi / 2.0;
		}
		else
		{
			return 1.5 * ec_fPi;
		}
	}
	else if((i_fImaginaryPart <= 1e-20) && (i_fImaginaryPart >= -1e-20) /*|| abs(i_fImaginaryPart) < 1e-10*/)
	{
		if(i_fRealPart >= 0.0)
		{
			return 0.0;
		}
		else
		{
			return ec_fPi;
		}
	}
	/*1st quadrant*/
	else if(i_fRealPart > 0.0 && i_fImaginaryPart > 0.0)
	{
		return atan(i_fImaginaryPart / i_fRealPart);
	}
	/*2nd*/
	else if(i_fRealPart < 0.0 && i_fImaginaryPart > 0.0)
	{
		return ec_fPi - atan(-i_fImaginaryPart / i_fRealPart);
	}
	/*3rd*/
	else if(i_fRealPart < 0.0 && i_fImaginaryPart < 0.0)
	{
		return ec_fPi + atan(i_fImaginaryPart / i_fRealPart);
	}
	/*4th*/
	else if(i_fRealPart > 0.0 && i_fImaginaryPart < 0.0)
	{
		return 2.0 * ec_fPi - atan(-i_fImaginaryPart / i_fRealPart);
	}
	else
	{
		std::cout << "OpticalElement::GetPhaseOfComplexNumber: Error calculating the phase of complex number " << i_fRealPart << " + " << i_fImaginaryPart << "i, set phase to 0.0" << std::endl;
		return 0.0;
	}
}

bool OpticalElement::PermissionToResetPointer(int i_nStackSize)
{
	if(m_nStackSizeBefore == i_nStackSize)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void OpticalElement::EndSimulation(int i_nNumberOfHistories)
{

}
