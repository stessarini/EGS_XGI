/*
###############################################################################
#
#   EGS_XGI Example_Usercode_Score_Energy
#   An example usercode for using EGS_XGI features within EGSnrc with energy
#		scoring features (compare EGSnrc tutor7pp)
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

#include "Example_Usercode_Score_Energy.h"

#include "xgi_global_variables.h"

#include "egs_interpolator.h"





Example_Usercode_Score_Energy::Example_Usercode_Score_Energy(int argc, char **argv)
	:EGS_XGIApplication(argc, argv),
	m_nNumberOfRayTracingPropagations(0),
	m_nNumberOfRayleighEvents(0)
{
	m_nNumberOfPhotonScatteringEvents = 0;

	m_nNumberOfRegions = geometry->regions();
	m_pDepositedEnergy = new EGS_ScoringArray(m_nNumberOfRegions + 3);
	m_fScoredEnergyThisHistory = 0.0;
	the_epcont->iausfl[ExtraEnergy] = 1;
}

Example_Usercode_Score_Energy::~Example_Usercode_Score_Energy()
{
	delete m_pDepositedEnergy;
}

int Example_Usercode_Score_Energy::ausgab(int iarg)
{
	int nStackAddress = the_stack->np-1;
	int nCharge = the_stack->iq[nStackAddress];
	bool bIsPrimary = e_bPrimary[nStackAddress];

	//score energy:
	if(iarg <= 4)
	{
		int nRegionNumber = the_stack->ir[nStackAddress] - 1;
		if(nRegionNumber == 0 && m_pSplittingAlgorithm->WhichOpticalElement() == "detector1")
		{
			nRegionNumber = m_nNumberOfRegions + 1;
		}
		EGS_Float aux = the_epcont->edep * the_stack->wt[nStackAddress];
		if(aux > 0.0)
		{
			m_pDepositedEnergy->score(nRegionNumber, aux);
			m_fScoredEnergyThisHistory += aux;
		}
	}

	if(iarg == UserDiscard || iarg == PegsCut || iarg == EgsCut)
	{
		//reset all the parameters. i.e. check if propagated all paths produced by a grating
		m_pSplittingAlgorithm->PotentialParameterReset();
		e_fPhase[nStackAddress] = 0.0;
		e_fLogNorm[nStackAddress] = 0.0;
		e_bPrimary[nStackAddress] = false;
		the_stack->latch[nStackAddress] = 0;
	}
	else if(iarg == BeforeTransport && nCharge == 0)
	{
		m_fPositionBeforeTVStep.x = the_stack->x[nStackAddress];
		m_fPositionBeforeTVStep.y = the_stack->y[nStackAddress];
		m_fPositionBeforeTVStep.z = the_stack->z[nStackAddress];
		//track phase  and norm
		int nCurrentRegion = the_stack->ir[nStackAddress] -2;
		if(bIsPrimary == true && nCurrentRegion > -1)
		{
			int nMediumIndex = geometry->medium(nCurrentRegion);
			EGS_Float fEnergy = the_stack->E[nStackAddress];
			EGS_Float steplength = the_epcont->tvstep;
			e_fPhase[nStackAddress] += steplength * (1.0 - m_pRefractiveIndexCalculator->GetRefractiveIndexDelta(nMediumIndex, fEnergy));
			if(nMediumIndex >= 0)
			{
				e_fLogNorm[nStackAddress] -= steplength * m_fMu_over_two_ForInitialEnergy[nMediumIndex];
			}
		}
	}
	else if(iarg == AfterTransport)
	{
		int nPreviousRegion =  the_epcont->irold - 2;
		int nCurrentRegion = the_stack->ir[nStackAddress] - 2;
		if(nCurrentRegion != nPreviousRegion)
		{
			m_pSplittingAlgorithm->PotentialParticleSplitting();
			if(nStackAddress == the_stack->np-1 && nCharge == 0 && nPreviousRegion > -1 && nCurrentRegion > -1)
			{
				int nMediumIndexNew = geometry->medium(nCurrentRegion);
				int nMediumIndexOld = geometry->medium(nPreviousRegion);
				if(nMediumIndexNew !=  nMediumIndexOld)
				{
					EGS_Vector fPosition;
					fPosition.x = the_stack->x[nStackAddress];
					fPosition.y = the_stack->y[nStackAddress];
					fPosition.z = the_stack->z[nStackAddress];

					EGS_Vector fDirection;
					fDirection.x = the_stack->u[nStackAddress];
					fDirection.y = the_stack->v[nStackAddress];
					fDirection.z = the_stack->w[nStackAddress];

					EGS_Float fEnergy = the_stack->E[nStackAddress];
					int nReflectionTransmission;
					nReflectionTransmission = ApplySnellsLawHowfar(&fPosition, &m_fPositionBeforeTVStep, &fDirection, nCurrentRegion, nPreviousRegion,  nMediumIndexNew, nMediumIndexOld, fEnergy);
					if(nReflectionTransmission == 0 || nReflectionTransmission == 1)
					{
						the_stack->u[nStackAddress] = fDirection.x;
						the_stack->v[nStackAddress] = fDirection.y;
						the_stack->w[nStackAddress] = fDirection.z;
					}
				}
			}
		}
	}
	else if(iarg == AfterPhoto || iarg == AfterCompton )
	{
		int nStackAddressOld = the_stack->npold-1;
		for(int nAddress = nStackAddressOld; nAddress <= nStackAddress; nAddress++)
		{
			e_fPhase[nAddress] = 0.0;
			e_fLogNorm[nAddress] = 0.0;
			e_bPrimary[nAddress] = false;
		}
	}
	else if(bIsPrimary == true)
	{
		if( iarg == BeforePhoto || iarg == BeforeCompton || iarg == BeforeRayleigh)
		{
			m_nInteractionCheck = iarg;
			int nStackSizeBeforeRayTracing = the_stack->np;

			EGS_Float fXPosition = the_stack->x[nStackAddress];
			EGS_Float fYPosition = the_stack->y[nStackAddress];
			EGS_Float fZPosition = the_stack->z[nStackAddress];

			EGS_Float fUDirection = the_stack->u[nStackAddress];
			EGS_Float fVDirection = the_stack->v[nStackAddress];
			EGS_Float fWDirection = the_stack->w[nStackAddress];

			int nRegionBeforeRayTracing = the_stack->ir[nStackAddress];
			int nRegionNewBerforeRayTracing = the_epcont->irnew;



			int nRTStackAddress = nStackAddress + 1;
			the_stack->x[nRTStackAddress] = the_stack->x[nStackAddress];
			the_stack->y[nRTStackAddress] = the_stack->y[nStackAddress];
			the_stack->z[nRTStackAddress] = the_stack->z[nStackAddress];
			the_stack->u[nRTStackAddress] = the_stack->u[nStackAddress];
			the_stack->v[nRTStackAddress] = the_stack->v[nStackAddress];
			the_stack->w[nRTStackAddress] = the_stack->w[nStackAddress];
			the_stack->ir[nRTStackAddress] = nRegionBeforeRayTracing;
			the_stack->E[nRTStackAddress] = the_stack->E[nStackAddress];
			the_stack->wt[nRTStackAddress] = the_stack->wt[nStackAddress];
			the_stack->latch[nRTStackAddress] = the_stack->latch[nStackAddress];
			e_fPhase[nRTStackAddress] = e_fPhase[nStackAddress];
			e_fLogNorm[nRTStackAddress] = e_fLogNorm[nStackAddress];
			e_bPrimary[nRTStackAddress] = true;
			the_stack->iq[nRTStackAddress] = 0;
			the_stack->np = nStackSizeBeforeRayTracing + 1;



			PerformRayTracing();
			//might have changed during Ray-Tracing
			the_epcont->irnew = nRegionNewBerforeRayTracing;
			the_stack->ir[nStackAddress] = nRegionBeforeRayTracing;

			if(the_stack->x[nStackAddress] != fXPosition || fUDirection != the_stack->u[nStackAddress])
			{
				std::cout << "warning RT changed the_stack" << std::endl;
			}
			if(iarg == BeforeRayleigh)
			{
				if(the_stack->latch[nStackAddress] == 0)
				{
					m_nNumberOfRayleighEvents++;
					e_bPrimary[nStackAddress] = true;
					the_stack->latch[nStackAddress] = BeforeRayleigh;
				}
				else
				{
					m_nNumberOfRayleighEvents++;
					e_fPhase[nStackAddress] = 0.0;
					e_fLogNorm[nStackAddress] = 0.0;
					e_bPrimary[nStackAddress] = false;
				}
			}
			else
			{
				//prepare for conventional MC transport
				m_nNumberOfPhotonScatteringEvents++;
				e_fPhase[nStackAddress] = 0.0;
				e_fLogNorm[nStackAddress] = 0.0;
				e_bPrimary[nStackAddress] = false;
			}
			if(nStackSizeBeforeRayTracing != the_stack->np)
			{
				std::cout << "Ray-tracing changed the_stack size by " << the_stack->np - nStackSizeBeforeRayTracing << std::endl;
				std::cout << "np: " << the_stack->np << std::endl;
				std::cout << "nStackSizeBeforeRayTracing: " << nStackSizeBeforeRayTracing << std::endl;
				std::cout << "Direction before ray-tracing: (" << fUDirection << ", " << fVDirection << ", " << fWDirection << ")" << std::endl;
				std::cout << "Position before ray-tracing: (" << fXPosition << ", " << fYPosition << ", " << fZPosition << ")" << std::endl;
			}
		}
	}
	return 0;
}

void Example_Usercode_Score_Energy::reportResults()
{
	egsInformation("\n\n last case = %d Etot = %g\n", (int)last_case,m_fEtot);
	std::cout << "Number of compton and photo effect scattering events of primary photons: " << m_nNumberOfPhotonScatteringEvents << std::endl;
	std::cout << "Number of Rayleigh scattering events of primary photons: " <<  m_nNumberOfRayleighEvents << std::endl;
	std::cout << "Pure ray tracing propagations: " << m_nNumberOfRayTracingPropagations << std::endl;
	std::cout << "Histories with too many reflections: " << m_nNumberOfParticlesWithToManyReflections <<std::endl;
	m_pSplittingAlgorithm->ReportSplittingSummary((int)last_case);

	std::cout << "ScoringArray output" << std::endl;
	EGS_Float norm = 1.0;
	m_pDepositedEnergy->reportResults(norm, "Deposited energy", false);
}


void Example_Usercode_Score_Energy::startHistory(EGS_I64 this_case)
{
	m_pDepositedEnergy->setHistory(this_case);
	if(m_fEnergyLastHistory != the_stack->E[0])
	{
		m_pSplittingAlgorithm->StartRayTracing(the_stack->E[0]);
		for(unsigned int i = 0; i < geometry->nMedia(); i++)
		{
			m_fMu_over_two_ForInitialEnergy[i] = 1.0 / i_gmfp[i].interpolate(log(the_stack->E[0])) / 2.0;
		}
	}
}

void Example_Usercode_Score_Energy::endHistory()
{
	m_pSplittingAlgorithm->endHistory();

	//make sure energy is conseved on history base and score energy in extra entry:
	m_pDepositedEnergy->score(m_nNumberOfRegions  + 2, abs(m_fScoredEnergyThisHistory - m_fEnergyLastHistory));

	m_fScoredEnergyThisHistory = 0.0;
	m_nInteractionCheck = 0;
}

int main(int argc, char **argv)
{
	Example_Usercode_Score_Energy simulation(argc, argv);

	simulation.run();

	simulation.reportResults();

	simulation.finish();

	return 0;
}
