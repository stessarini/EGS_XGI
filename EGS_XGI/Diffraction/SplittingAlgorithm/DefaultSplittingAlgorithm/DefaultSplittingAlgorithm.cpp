/*
###############################################################################
#
#   EGS_XGI DefaultSplittingAlgorithm
#   Default interferometer class, managing splitting of the optics components
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
#include "DefaultSplittingAlgorithm.h"



#include "DetectorZPlane.h"

#include "xgi_global_variables.h"


DefaultSplittingAlgorithm::DefaultSplittingAlgorithm()
        :SplittingAlgorithm()
{
}


DefaultSplittingAlgorithm::~DefaultSplittingAlgorithm()
{

}

void DefaultSplittingAlgorithm::ReportSplittingAlgorithm()
{
	std::cout << "+++++++++++++++++++++++++++++++" << endl;
	std::cout << "+++SplittingAlgorithm Report+++" << endl;
	std::cout << "Type: DefaultSplittingAlgorithm"<<endl;
	std::cout << "source coherence: ";
	if(m_bCoherentSource)
	{
		std::cout << "coherent" << std::endl;
	}
	else
	{
		std::cout << "incoherent" << std::endl;
	}
	std::cout << "Number of OpticalElements: " << m_pSplittingObjects.size() << endl << endl;
	std::cout << "---Included OpticalElements---" << endl;
	for(vector<OpticalElement*>::iterator SOIterator = m_pSplittingObjects.begin();  SOIterator != m_pSplittingObjects.end();  SOIterator++)
	{
		(*SOIterator)->ReportOpticalElement();
		std::cout << "-------------------------------" << endl;
	}

}


void DefaultSplittingAlgorithm::ReportSplittingSummary(int ncase)
{
	unsigned int nTotalNumberOfGeneratedHistories;
	vector<unsigned int> SplittingEvents;
	vector<unsigned int> GeneratedClones;
	std::cout << "+++++++++++++++++++++++++++++++" << endl;
	std::cout << "+++SplittingAlgorithm Summary++" << endl;
	std::cout << "OpticalElements:" << endl;
	for(vector<OpticalElement*>::iterator SOIterator = m_pSplittingObjects.begin();  SOIterator != m_pSplittingObjects.end();  SOIterator++)
	{
		(*SOIterator)->PrintSummary();
		std::cout << "-------------------------------" << endl;
		SplittingEvents.push_back((*SOIterator)->GetNumberOfTimesAppliedOpticalElemetRule());
	}
	std::cout << "-------------------------------" << endl;
	std::cout << "Started histories = " << ncase << endl;
	std::cout << "+++++++++++++++++++++++++++++++" << endl;
}


//Assume that it was checked that its a photon && irnew!=-1 && it is called after transport
void DefaultSplittingAlgorithm::PotentialParticleSplitting()
{

	if(m_pLatestActiveSplittingObject != m_pSplittingObjects.begin())
	{
		vector<OpticalElement*>::iterator NextSplittingObject = m_pLatestActiveSplittingObject - 1;
		if((*NextSplittingObject)->ApplyOpticalElementRule())//Splitting was applied
		{
			//Update the pointer
			m_pLatestActiveSplittingObject = NextSplittingObject;
		}
	}
}

//Reset Parameters if necessary
//Assume that it is given that: ausgab::iarg = After_discard
void DefaultSplittingAlgorithm::PotentialParameterReset()
{
	if(!(m_pLatestActiveSplittingObject == m_pSplittingObjects.end()))
	{
		//Check whether there are stil some split paths from m_pLatestActiveSplittingObject to be processed
    while((m_pLatestActiveSplittingObject != m_pSplittingObjects.end()) && (*m_pLatestActiveSplittingObject)->PermissionToResetPointer(the_stack->np -1) )
    {
    	m_pLatestActiveSplittingObject++;//Set the pointer to the previous OpticalElement
    }
	}
	//else nothing to do
}


int DefaultSplittingAlgorithm::InitDerivedSplittingAlgorithm(EGS_Input* i_Input)
{
	std::vector<std::string> sSetup;
	std::vector<std::string> sSourceCoherence;

	int errSetup = i_Input->getInput("setup", sSetup);
	int errCoherence = i_Input->getInput("source coherence", sSourceCoherence);
	//In case coherence is not set in the input file a warning should be shown. In that case an incoherent source is assumed
	if(errCoherence)
	{
		std::cout << "DefaultSplittingAlgorithm::InitDerivedSplittingAlgorithm: Error reading 'coherence'-key." << std::endl;
		std::cout << "DefaultSplittingAlgorithm::InitDerivedSplittingAlgorithm: Set to 'incoherent'." << std::endl;
		m_bCoherentSource = false;
	}
	else
	{
		if(sSourceCoherence[0] == "coherent")
		{
			m_bCoherentSource = true;
		}
		else if(sSourceCoherence[0] == "incoherent")
		{
			m_bCoherentSource = false;
		}
		else
		{
			std::cout << "DefaultSplittingAlgorithm::InitDerivedSplittingAlgorithm: Unknown 'coherence'-key value" << std::endl;
			m_bCoherentSource = false;
		}
	}

	if(errSetup)
	{
		std::cout << "DefaultSplittingAlgorithm::InitDerivedSplittingAlgorithm: Error reading 'setup'." << std::endl;
		return 0;
	}
	else
	{
		if(sSetup.size()>0)
		{
			if(m_pAllSplittingObjectsInInputfile.size() != sSetup.size())
			{
				std::cout << "DefaultSplittingAlgorithm::InitDerivedSplittingAlgorithm: Warning setup doesn't contain all optical elements defined in the input file..." << std::endl;
			}
			std::vector<OpticalElement*> Setup;
			for(int nSetupCounter = sSetup.size()-1; nSetupCounter != -1; nSetupCounter--)
			{
				OpticalElement* NextOpticalElement = GetPointerTo(sSetup[nSetupCounter]);
				if(NextOpticalElement != NULL)
				{
					Setup.push_back(NextOpticalElement);
				}
				else
				{
					std::cout << "DefaultSplittingAlgorithm::InitDerivedSplittingAlgorithm: Error optical element with name " << sSetup[nSetupCounter] << "doesn't exist" << std::endl;
					return 0;
				}
			}
			m_pSplittingObjects = Setup;
			m_pLatestActiveSplittingObject = m_pSplittingObjects.end();
		}
		else
		{
			std::cout << "DefaultSplittingAlgorithm::InitDerivedSplittingAlgorithm: No setup defined." << std::endl;
			return 0;
		}
	}
  return 1;
}


void DefaultSplittingAlgorithm::endHistory()
{
	if(m_bCoherentSource == false)
	{
		for (vector<OpticalElement*>::iterator OpticalElementIterator = m_pSplittingObjects.begin(); OpticalElementIterator != m_pSplittingObjects.end(); OpticalElementIterator++)
		{
			(*OpticalElementIterator)->EndRayTracing();
		}
	}
}


void DefaultSplittingAlgorithm::EndSimulation(int i_nNumberOfHistories)
{
  if(m_bCoherentSource == true)
  {
    for (vector<OpticalElement*>::iterator OpticalElementIterator = m_pSplittingObjects.begin(); OpticalElementIterator != m_pSplittingObjects.end(); OpticalElementIterator++)
		{
			(*OpticalElementIterator)->EndRayTracing();
		}
  }
  for (vector<OpticalElement*>::iterator OpticalElementIterator = m_pSplittingObjects.begin(); OpticalElementIterator != m_pSplittingObjects.end(); OpticalElementIterator++)
  {
    (*OpticalElementIterator)->EndSimulation(i_nNumberOfHistories);
  }
}


std::string DefaultSplittingAlgorithm::WhichOpticalElement()
{
  if(m_pLatestActiveSplittingObject != m_pSplittingObjects.end())
  {
    return (*m_pLatestActiveSplittingObject)->GetName();
  }
  else
  {
    return "DefaultSplittingAlgorithm::source";
  }
}


int DefaultSplittingAlgorithm::IsInInitialCondition()
{
  if(m_pLatestActiveSplittingObject != m_pSplittingObjects.end())
  {
    //something went wrong during splitting
    return 1;
  }
  else
  {
    return 0;
  }
}
