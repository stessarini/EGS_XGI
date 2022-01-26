/*
###############################################################################
#
#   EGS_XGI SplittingObject
#   Base class for splitting objects such as gratings
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
#include "SplittingObject.h"


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
