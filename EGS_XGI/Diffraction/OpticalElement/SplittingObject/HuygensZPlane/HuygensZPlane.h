/*
###############################################################################
#
#   EGS_XGI HuygensZPlane header
#   Apply uniform path splitting on a plane orthogonal to z-axis with anglular
#   restrictions.
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
#ifndef _HUYGENSZPLANE_H_
#define _HUYGENSZPLANE_H_

#include "SplittingObject.h"
#include "egs_rndm.h"

class HuygensZPlane : public SplittingObject
{
public:
  HuygensZPlane(std::string i_sName, EGS_RandomGenerator* i_pRandomNumberGenerator);
  virtual ~HuygensZPlane();

  virtual bool ApplyOpticalElementRule();
	virtual void PrintSummary();
	virtual void ReportOpticalElement();
  virtual int InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry);

protected:

  EGS_Float m_fZPosition;
  EGS_Float m_fXmin, m_fXmax;
  EGS_Float m_fYmin, m_fYmax;

  //minimum and maximum angle of split paths w.r.t optical (z) axis
  EGS_Float m_fSplittingAngleMin, m_fSplittingAngleMax;
  EGS_Float m_fRangePolarAngle;

  EGS_RandomGenerator* m_pRandomNumberGenerator;


};
#endif
