/*
###############################################################################
#
#   EGS_XGI HuygensZPlane
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
#include "HuygensZPlane.h"
#include "xgi_global_variables.h"

#include "egs_interface2.h"

HuygensZPlane::HuygensZPlane(std::string i_sName, EGS_RandomGenerator* i_pRandomNumberGenerator)
  :SplittingObject(i_sName, "HuygensZPlane"),
  m_fZPosition(0.0),
  m_fXmin(0.0),
  m_fXmax(0.0),
  m_fYmin(0.0),
  m_fYmax(0.0),
  m_fSplittingAngleMin(0.0),
  m_fSplittingAngleMax(0.0),
  m_pRandomNumberGenerator(i_pRandomNumberGenerator)
{
  m_nSplittingNumber = 1;
}

///////////////////////////////////


HuygensZPlane::~HuygensZPlane()
{

}

///////////////////////////////////


bool HuygensZPlane::ApplyOpticalElementRule()
{
  m_nStackSizeBefore = the_stack->np-1;
  EGS_Float fZPosition = the_stack->z[m_nStackSizeBefore];
  if(abs(m_fZPosition - fZPosition) < 1e-10)
  {
    EGS_Float fXPosition = the_stack->x[m_nStackSizeBefore];
    if(m_fXmin < fXPosition && fXPosition < m_fXmax)
    {
      EGS_Float fYPosition = the_stack->y[m_nStackSizeBefore];
      if(m_fYmin < fYPosition && fYPosition < m_fYmax)
      {
        EGS_I32 nRegionIndexBefore = the_stack->ir[m_nStackSizeBefore];
        if(e_bPrimary[m_nStackSizeBefore] == true && nRegionIndexBefore > 1 && the_stack->iq[m_nStackSizeBefore] == 0)
        {
          //Get all parameters before the splitting event:
          EGS_Float fWeightAfterPathSplitting = the_stack->wt[m_nStackSizeBefore] / m_nSplittingNumber;
  				EGS_Float fPhaseBeforePathSplitting = e_fPhase[m_nStackSizeBefore];
  				EGS_Float fNormBeforeSpitting = e_fLogNorm[m_nStackSizeBefore];

  				EGS_Float fEnergyBefore = the_stack->E[m_nStackSizeBefore];

          EGS_Vector fDirectionBefore;
          fDirectionBefore.x= the_stack->u[m_nStackSizeBefore];
          fDirectionBefore.y= the_stack->v[m_nStackSizeBefore];
          fDirectionBefore.z= the_stack->w[m_nStackSizeBefore];
          fDirectionBefore.normalize();
          EGS_Float fRadiusPolar = sqrt(1.0 - fDirectionBefore.y * fDirectionBefore.y);

  				EGS_Float dnearBefore = the_stack->dnear[m_nStackSizeBefore];
  				EGS_I32 nLatchBerfore = the_stack->latch[m_nStackSizeBefore];

          for(int nPathCounter = 0; nPathCounter < m_nSplittingNumber; nPathCounter++)
          {
            //generate new random direction
            //Assuming grating lines along y-direction, i.e., periodicity along x-direction
            //Get polar angle of new direction vector
            EGS_Float fPolarAngle = m_pRandomNumberGenerator->getUniform() * m_fRangePolarAngle + m_fSplittingAngleMin;

            e_fPhase[m_nStackSizeBefore + nPathCounter] = fPhaseBeforePathSplitting;
  					e_fLogNorm[m_nStackSizeBefore + nPathCounter] = fNormBeforeSpitting;
  					e_bPrimary[m_nStackSizeBefore + nPathCounter] = true;

  					the_stack->E[m_nStackSizeBefore + nPathCounter] 	= fEnergyBefore;

  					the_stack->x[m_nStackSizeBefore + nPathCounter] 	= fXPosition;
  					the_stack->y[m_nStackSizeBefore + nPathCounter] 	= fYPosition;
  					the_stack->z[m_nStackSizeBefore + nPathCounter] 	= fZPosition;

  					the_stack->u[m_nStackSizeBefore + nPathCounter] 	= fRadiusPolar * sin(fPolarAngle);
  					the_stack->v[m_nStackSizeBefore + nPathCounter] 	= fDirectionBefore.y;
  					the_stack->w[m_nStackSizeBefore + nPathCounter] 	= fRadiusPolar * cos(fPolarAngle);

  					the_stack->dnear[m_nStackSizeBefore + nPathCounter] = dnearBefore;

  					the_stack->wt[m_nStackSizeBefore + nPathCounter] 	= fWeightAfterPathSplitting;

  					the_stack->iq[m_nStackSizeBefore + nPathCounter] 	= (EGS_I32)0;
  					the_stack->ir[m_nStackSizeBefore + nPathCounter] 	= nRegionIndexBefore;
  					the_stack->latch[m_nStackSizeBefore + nPathCounter] = nLatchBerfore;
          }
          the_stack->np += m_nSplittingNumber -1;
          m_nNumberOfGeneratedPaths += m_nSplittingNumber -1;
  				m_nNumberOfTimesApplied++;
          return true;
        }
        else
        {
          //only falg that it went through *this
          return true;
        }
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}

///////////////////////////////////


void HuygensZPlane::PrintSummary()
{
  std::cout << "--------------------------------" << std::endl;
	std::cout << "Splitting Object of Type: " << m_sType << std::endl;
	std::cout << "Name: " << m_sName << std::endl;
	std::cout << "Number of splittings: " << m_nNumberOfTimesApplied << std::endl;
	std::cout << "Number of paths created: " << m_nNumberOfGeneratedPaths << std::endl;
	std::cout << "--------------------------------" << std::endl;
}

///////////////////////////////////


void HuygensZPlane::ReportOpticalElement()
{
  std::cout << "--------------------------------" << std::endl;
  std::cout << "Splitting Object of Type: HuygensZPlane" << std::endl;
  std::cout << "Name: " << m_sName << std::endl;
  std::cout << "Splitting number = " << m_nSplittingNumber << std::endl;
  std::cout << "z-position: " << m_fZPosition << std::endl;
  std::cout << "[Xmin, Xmax] = [" << m_fXmin << ", " << m_fXmax << "]" << std::endl;
  std::cout << "[Ymin, Ymax] = [" << m_fYmin << ", " << m_fYmax << "]" << std::endl;
  std::cout << "angle range : " << m_fRangePolarAngle << std::endl;
  std::cout << "angle min : " << m_fSplittingAngleMin << std::endl;
  std::cout << "--------------------------------" << std::endl;
}

///////////////////////////////////


int HuygensZPlane::InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry)
{
  EGS_Float fPosition;
  std::vector<EGS_Float> fXRange;
  std::vector<EGS_Float> fYRange;
  std::vector<EGS_Float> fPolarRange;
  int nSplittingNumber;

  int errPos = i_sInput->getInput("position", fPosition);
  int errXR = i_sInput->getInput("x-limits", fXRange);
  int errYR = i_sInput->getInput("y-limits", fYRange);
  int errPhiR = i_sInput->getInput("angle limits", fPolarRange);
  int errSN = i_sInput->getInput("number of paths", nSplittingNumber);

  if(errPos)
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: " << std::endl;
    return 0;
  }
  else
  {
    m_fZPosition = fPosition;
  }

  if(errXR)
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: reading 'x-limits'." << std::endl;
    return 0;
  }
  else if(fXRange.size() != 2)
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: unexpected format of 'x-limits'. Requires exactly 2 inputs: Xmin and Xmax." << std::endl;
    return 0;
  }
  else if(fXRange[0]>=fXRange[1])
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: 'x-limits' has to be given in increasing order." << std::endl;
    return 0;
  }
  else
  {
    m_fXmin = fXRange[0];
    m_fXmax = fXRange[1];
  }

  if(errYR)
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: reading 'y-limits'." << std::endl;
    return 0;
  }
  else if(fYRange.size() != 2)
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: unexpected format of 'y-limits'. Requires exactly 2 inputs: Ymin and Ymax." << std::endl;
    return 0;
  }
  else if(fYRange[0]>=fYRange[1])
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: 'y-limits' has to be given in increasing order." << std::endl;
    return 0;
  }
  else
  {
    m_fYmin = fYRange[0];
    m_fYmax = fYRange[1];
  }

  if(errPhiR)
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: reading 'angle limits'." << std::endl;
    return 0;
  }
  else if(fPolarRange.size() != 2)
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: unexpected format of ' angle limits'. Requires exactly 2 inputs: min angle and max angle." << std::endl;
    return 0;
  }
  else if(fPolarRange[0]>=fPolarRange[1])
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: 'angle limits' has to be given in increasing order." << std::endl;
    return 0;
  }
  else
  {
    m_fSplittingAngleMin = fPolarRange[0];
    m_fSplittingAngleMax = fPolarRange[1];
    m_fRangePolarAngle = m_fSplittingAngleMax - m_fSplittingAngleMin;
    if(m_fSplittingAngleMax >= ec_fPi / 2.0 || m_fSplittingAngleMin <= -ec_fPi / 2.0)
    {
      std::cout << "HuygensZPlane::InitOpticalElement: Warning: 'angle limits' outside (-pi/2, pi/2)." << std::endl;
    }
  }

  if(errSN)
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: reading 'number of paths'." << std::endl;
    return 0;
  }
  if(nSplittingNumber <= 0)
  {
    std::cout << "HuygensZPlane::InitOpticalElement: Error: 'number of paths' has to be greater than 0." << std::endl;
    return 0;
  }
  else
  {
    m_nSplittingNumber = nSplittingNumber;
  }

  return 1;

}

///////////////////////////////////
