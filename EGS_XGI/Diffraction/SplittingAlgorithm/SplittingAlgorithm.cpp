#include "SplittingAlgorithm.h"

#include "DetectorZPlane.h"
#include "DetectorZPlaneMem.h"

#include "FourierSeriesZPlane.h"
#include "BinarySourceGratingFourierZPlane.h"
#include "HuygensZPlane.h"

#include "egs_interpolator.h"

SplittingAlgorithm::SplittingAlgorithm()
{
}

SplittingAlgorithm::~SplittingAlgorithm()
{
	for(unsigned int nOpticalElementCounter = 0; nOpticalElementCounter < m_pAllSplittingObjectsInInputfile.size(); nOpticalElementCounter++)
	{
		if(m_pAllSplittingObjectsInInputfile[nOpticalElementCounter] != NULL)
		{
			delete m_pAllSplittingObjectsInInputfile[nOpticalElementCounter];
			m_pAllSplittingObjectsInInputfile[nOpticalElementCounter] = NULL;
		}
	}
}

SplittingAlgorithm::SplittingAlgorithm(std::string i_sInputFile)
{
}


void SplittingAlgorithm::ReportSplittingAlgorithm()
{
	std::cout << "----------------------------------" << endl;
	std::cout << "Report on the splitting algorithm:" << endl;
}

void SplittingAlgorithm::ReportSplittingSummary(int ncase)
{
	std::cout<<"Report on the splitting algorithm after run:" << endl;
}


void SplittingAlgorithm::PotentialParticleSplitting()
{
}

void SplittingAlgorithm::PotentialParameterReset()
{
}

void SplittingAlgorithm::EndSimulation(int i_nNumberOfHistories)
{
}

int SplittingAlgorithm::InitSplittingAlgorithm(EGS_Input* i_Input, std::string* i_sInputFileName, EGS_RandomGenerator* i_RandomNumberGenerator, EGS_BaseGeometry* i_pGeometry, EGS_Interpolator* i_gmfp, RefractiveIndexCalculator* i_pRefractiveIndexCalculator)
{
	//Import the parameters for the OpticalElements
	EGS_Input* OpticalElementInput;
	OpticalElementInput = i_Input->takeInputItem("optical element");

	while(OpticalElementInput != 0)
	{
		std::vector<std::string> sType;
		std::vector<std::string> sName;
		int errType = OpticalElementInput->getInput("type", sType);
		int errName = OpticalElementInput->getInput("name", sName);
		if(errType || errName)
		{
			std::cout << "SplittingAlgorithm::InitSplittingAlgorithm: Error: no type or name specified."  << std::endl;
			return 0;
		}
		else if(sType[0] == "DetectorZPlane")
		{
			OpticalElement* pNewElement = new DetectorZPlane(sName[0], *i_sInputFileName);
			m_pAllSplittingObjectsInInputfile.push_back(pNewElement);
			int errNewElement = pNewElement->InitOpticalElement(OpticalElementInput, i_pGeometry);
			if(errNewElement == 0)
			{
				return 0;
			}
		}
		else if (sType[0] == "FourierSeriesZPlane")
		{
			OpticalElement* pNewElement = new FourierSeriesZPlane(sName[0], i_gmfp, i_pRefractiveIndexCalculator);
			m_pAllSplittingObjectsInInputfile.push_back(pNewElement);
			int errNewElement = pNewElement->InitOpticalElement(OpticalElementInput, i_pGeometry);
			if (errNewElement == 0)
			{
				return 0;
			}
		}
		else if(sType[0] == "DetectorZPlaneMem")
		{
			OpticalElement* pNewElement = new DetectorZPlaneMem(sName[0], *i_sInputFileName);
			m_pAllSplittingObjectsInInputfile.push_back(pNewElement);
			int errNewElement = pNewElement->InitOpticalElement(OpticalElementInput, i_pGeometry);
			if(errNewElement == 0)
			{
				std::cout << "SplittingAlgorithm: DetectorZPlaneMem" << std::endl;
				return 0;
			}
		}
		else if(sType[0] == "BinarySourceGratingFourierZPlane")
		{
			OpticalElement* pNewElement = new BinarySourceGratingFourierZPlane(sName[0], i_RandomNumberGenerator);
			m_pAllSplittingObjectsInInputfile.push_back(pNewElement);
			int errNewElement = pNewElement->InitOpticalElement(OpticalElementInput, i_pGeometry);
			if(errNewElement == 0)
			{
				std::cout << "SplittingAlgorithm: Error setting up BinarySourceGratingFourierZPlane" << std::endl;
				return 0;
			}
		}
		else if(sType[0] == "HuygensZPlane")
		{
			OpticalElement* pNewElement = new HuygensZPlane(sName[0], i_RandomNumberGenerator);
			m_pAllSplittingObjectsInInputfile.push_back(pNewElement);
			int errNewElement = pNewElement->InitOpticalElement(OpticalElementInput, i_pGeometry);
			if(errNewElement == 0)
			{
				std::cout << "SplittingAlgorithm: Error setting up HuygensZPlane" << std::endl;
				return 0;
			}
		}
		else
		{
			std::cout << "Unknown optical element type.." << std::endl;
			return 0;
		}
		delete OpticalElementInput;
		OpticalElementInput = i_Input->takeInputItem("optical element");
	}
	return 1;
}

int SplittingAlgorithm::InitDerivedSplittingAlgorithm(EGS_Input* i_Input)
{
	return 0;
}

OpticalElement* SplittingAlgorithm::GetPointerTo(std::string i_sName)
{
	for(unsigned int OpticalElementCounter = 0; OpticalElementCounter != m_pAllSplittingObjectsInInputfile.size(); OpticalElementCounter++)
	{
		if(m_pAllSplittingObjectsInInputfile[OpticalElementCounter]->GetName() == i_sName)
		{
			return m_pAllSplittingObjectsInInputfile[OpticalElementCounter];
		}
	}
	return NULL;
}

void SplittingAlgorithm::endHistory()
{

}

void SplittingAlgorithm::StartRayTracing(EGS_Float i_fInitialEnergy)
{

}

std::string SplittingAlgorithm::WhichOpticalElement()
{
	return (* (m_pSplittingObjects.begin()))->GetName();
}


int SplittingAlgorithm::IsInInitialCondition()
{
	return 0;
}
