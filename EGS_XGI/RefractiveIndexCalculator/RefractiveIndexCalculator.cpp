/*
###############################################################################
#
#   EGS_XGI RefractiveIndexCalculator
#   Class to handle the decrement \delta of the real part of the refractive
#   index: n = 1 - \delta
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

#include "RefractiveIndexCalculator.h"
#include <map>
#include "xgi_global_variables.h"

/**********************/


RefractiveIndexCalculator::RefractiveIndexCalculator()
	:m_fRefractiveIndexData()
{

}


/**********************/



/**********************/


RefractiveIndexCalculator::~RefractiveIndexCalculator()
{

}


/**********************/


/*Assume that the real part of the refractive index is given by: $n = 1 + \delta$,
	with a quadratic dependence on the wavelength $\lambda$:
	$\delta = \delta^{*} * \lambda^2$.
	!This class assumes that the input file for the refractive index (.refind) has the following format
	<medium_name_1> <\delta_1^{*}>
	<medium_name_2> <\delta_2^{*}>
	...
	with the medium names matching the PEGS file.
	Note: Not all media defined in the PEGS file need to be listed. Only those used in the simulation run
		  are needed
	Note: $\delta^{*} is assumed to be 0 for all medium names that are not found in the refractive index file$*/
void RefractiveIndexCalculator::ImportRefractiveIndex(string i_sRefFile, EGS_BaseGeometry* i_pSimulationGeometry)
{
	fstream file;
	file.open(i_sRefFile);
	string sLine;

	unsigned int NumberOfMedia = i_pSimulationGeometry->nMedia();
	vector<string> sMediaNames;
	for (unsigned int i = 0; i < NumberOfMedia; i++)
	{
		sMediaNames.push_back(i_pSimulationGeometry->getMediumName(i));
	}
	//A list to store media for which there is no refractive index data.
	map<int, EGS_Float> Buffer;

	while (getline(file, sLine) && m_fRefractiveIndexData.size() < NumberOfMedia)
	{
		unsigned int nSpace = sLine.find(" ");
		if (nSpace != string::npos)
		{
			string copy = sLine;
			sLine.erase(nSpace, string::npos);

			for (unsigned int i = 0; i < sMediaNames.size(); i++)
			{
				if ((sMediaNames[i]).compare(sLine) == 0)
				{
					copy.erase(0, nSpace + 1);
					EGS_Float value = stod(copy);
					Buffer.insert(pair<int, EGS_Float>(i, value));
					break;
				}
			}
		}
	}

	for (unsigned int i = 0; i < NumberOfMedia; i++)
	{
		if (Buffer.find(i) != Buffer.end())
		{
			m_fRefractiveIndexData.push_back(Buffer[i]);
		}
		else
		{
			cout << "RefractiveIndexCalculator::ImportRefractiveIndex::WARNIG did not find the medium " << sMediaNames[i] << " in the .refind file." << endl;
			m_fRefractiveIndexData.push_back(0.0);
		}
	}

	file.close();
	for (unsigned int i = 0; i < NumberOfMedia; i++)
	{
		cout << "RefractiveIndexCalculator::m_fRefractiveIndexData: [" << i_pSimulationGeometry->getMediumName(i) << "] = " << m_fRefractiveIndexData[i] << endl;
	}

}


/**********************/



/**********************/


EGS_Float RefractiveIndexCalculator::GetRefractiveIndexDelta(int i_nMediumIndex, EGS_Float i_fEnergy)
{
	if (i_nMediumIndex >= 0)
	{
		return m_fRefractiveIndexData[i_nMediumIndex] * pow(ec_fEnergyToWaveLength / i_fEnergy, 2);
	}
	else
	{
		return 0.0;//VACCUM has a medium index of -1
	}
}


/**********************/



/**********************/


EGS_Float RefractiveIndexCalculator::GetRefractiveIndexDeltaFromWavelength(int i_nMediumIndex, EGS_Float i_fWavelength)
{
	if (i_nMediumIndex >= 0)
	{
		return m_fRefractiveIndexData[i_nMediumIndex] * i_fWavelength * i_fWavelength;
	}
	else
	{
		return 0.0;//VACCUM has a medium index of -1
	}
}


/**********************/
