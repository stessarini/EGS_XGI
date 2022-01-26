#include "SplittingObject.h"


//using namespace std;



/***********************************************************/
/* SplittingObject */
/***********************************************************/
SplittingObject::SplittingObject(std::string i_sName, std::string i_sType)
	:OpticalElement(),
	m_nMinimumStackSizeAfter(-1),
	m_sType(i_sType),
	m_nNumberOfGeneratedPaths(0)
{
	m_sName = i_sName;
}

SplittingObject::~SplittingObject()
{

}

string SplittingObject::GetName()
{
	return m_sName;
}


std::string SplittingObject::GetType()
{
	return m_sType;
}

int SplittingObject::GetNumberOfSplittingEvents()
{
	return m_nNumberOfTimesApplied;
}

int SplittingObject::GetNumberOfGeneratedPaths()
{
	return m_nNumberOfGeneratedPaths;
}

bool SplittingObject::IsOnOpticalElement (EGS_Vector* i_pPosition)
{
	return false;
}

bool SplittingObject::ApplyOpticalElementRule()
{
	return true;
}

bool SplittingObject::PermissionToResetPointer(int i_nStackSize)
{
	if(m_nStackSizeBefore == i_nStackSize)
	{
		return true;
	}
	else if(m_nStackSizeBefore >= i_nStackSize)
	{
		std::cout << "SplittingObject::PermissionToResetPointer: "<< m_sName << ": warning stack size too low: " << m_nStackSizeBefore << " > " << i_nStackSize << std::endl;
		return true;
	}
	else
	{
		return false;
	}
}

void SplittingObject::PrintSummary()
{
	std::cout << "--------------------------------" << std::endl;
	std::cout << "Splitting Object of Type: " << m_sType << std::endl;
	std::cout << "Name: " << m_sName << std::endl;
	std::cout << "Number of splittings: " << m_nNumberOfTimesApplied << std::endl;
	std::cout << "Number of paths created: " << m_nNumberOfGeneratedPaths << std::endl;
	std::cout << "--------------------------------" << std::endl;
}

void SplittingObject::ReportOpticalElement()
{
	std::cout << "Splitting Object of Type: " << m_sType << std::endl;
	std::cout << "Name: " << m_sName << std::endl;
	std::cout << "Splitting number = " << m_nSplittingNumber << std::endl;
}
