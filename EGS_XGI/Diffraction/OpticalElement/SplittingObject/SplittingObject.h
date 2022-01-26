/*
###############################################################################
#
#   EGS_XGI SplittingObject header
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
#ifndef _SPLITTINGOBJECT_H_
#define _SPLITTINGOBJECT_H_


#include <string>

#include "OpticalElement.h"


class SplittingObject : public OpticalElement
{
public:


	SplittingObject(std::string i_sName, std::string i_sType);
	virtual ~SplittingObject();

	virtual bool IsOnOpticalElement (EGS_Vector* i_pPosition);
	virtual bool ApplyOpticalElementRule();
	virtual void PrintSummary();
	virtual void ReportOpticalElement();

	virtual bool PermissionToResetPointer(int i_nStackSize);

	std::string GetName();
	std::string GetType();

	int GetNumberOfSplittingEvents();
	int GetNumberOfGeneratedPaths();


protected:

	int m_nMinimumStackSizeAfter;
  std::string m_sType;
  unsigned int m_nNumberOfGeneratedPaths;
  int m_nSplittingNumber;
};

#endif
