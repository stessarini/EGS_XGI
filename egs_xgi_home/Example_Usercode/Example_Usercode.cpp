/***************************************************
*	Example_Usercode.cpp
*		- Derive from EGS_AdvancedApplication
*
*
*
*
****************************************************/

#include "Example_Usercode.h"

#include "xgi_global_variables.h"

#include "egs_interpolator.h"





Example_Usercode::Example_Usercode(int argc, char **argv)
	:EGS_XGIApplication(argc, argv),
	m_nNumberOfRayTracingPropagations(0),
	m_nNumberOfRayleighEvents(0)
{
	m_nNumberOfPhotonScatteringEvents = 0;

	/************************************************************/
	//extract data for python:
	/*std::cout.precision(17);
	std::cout << "mu and delta at 20 keV of" << std::endl;
	for(int nMediaIndex = 0; nMediaIndex != geometry->nMedia(); nMediaIndex++)
	{
		std::cout << geometry->getMediumName(nMediaIndex) << ": " << 1.0 / i_gmfp[nMediaIndex].interpolate(log(0.020)) << std::endl;
		std::cout << "\t" << m_pRefractiveIndexCalculator->GetRefractiveIndexDelta(nMediaIndex, 0.020) <<std::endl;
	}
*/
	/*std::vector<EGS_Float> energies{0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041};
	std::cout << "mu of " << std::endl;
	for(int nMediaIndex = 0; nMediaIndex != geometry->nMedia(); nMediaIndex++)
	{
		std::cout << geometry->getMediumName(nMediaIndex) << ": ";
		for(int j = 0; j < energies.size(); j++)
		{
			std::cout << 1.0 / i_gmfp[nMediaIndex].interpolate(log(energies[j]))<<", ";
		}
		std::cout << std::endl;
	}*/
	/************************************************************/
}

Example_Usercode::~Example_Usercode()
{

}

int Example_Usercode::ausgab(int iarg)
{
	int nStackAddress = the_stack->np-1;
	int nCharge = the_stack->iq[nStackAddress];
	bool bIsPrimary = e_bPrimary[nStackAddress];
	//the_stack->latch[nStackAddress] = iarg;
//simple way of following setp by step:
/*if(bIsPrimary == false && nCharge == 0)
{*/
	/*std::cout << "------------------------------------------" << std::endl;
	std::cout << "iarg = " << iarg << std::endl;
	std::cout << "charge = " << nCharge << std::endl;
	std::cout << "......" << std::endl;
	std::cout << "np-1 = " << the_stack->np-1 << std::endl;
	std::cout << "the_stack->npold-1 = " << the_stack->npold - 1 << std::endl;
	std::cout << "......" << std::endl;
	std::cout << "r = (" << the_stack->x[nStackAddress] << ", " << the_stack->y[nStackAddress]  << ", " << the_stack->z[nStackAddress] <<")" << std::endl;
	std::cout << "p = " << the_stack->u[nStackAddress] << ", " << the_stack->v[nStackAddress] << ", " << the_stack->w[nStackAddress] << ")" << std::endl;
	std::cout << "e_fLogNorm[np-1] = " <<  "= " << e_fLogNorm[nStackAddress] << std::endl;
	std::cout << "......" << std::endl;
	std::cout << "the_stack->ir[nStackAddress] = " << the_stack->ir[nStackAddress] - 2 << std::endl;
	std::cout << "the_epcont->irold = " << the_epcont->irold - 2 << std::endl;
	std::cout << "the_epcont->irnew = " << the_epcont->irnew - 2 << std::endl;
	std::cout << "......" << std::endl;
	std::cout << "the_stack->wt[nStackAddress] = " << the_stack->wt[nStackAddress] << std::endl;
	std::cout << "e_fPhase[np-1] = " << e_fPhase[nStackAddress] << std::endl;
	std::cout << "inside geometry: " << geometry->getName() << std::endl;
	int imed = geometry->medium(the_stack->ir[nStackAddress] - 2);
	std::cout << "medium index: " << imed << std::endl;
	if(imed >= 0)
	{
		std::cout << "medium name: " << geometry->getMediumName(imed) << std::endl;
	}
	std::cout << m_pSplittingAlgorithm->WhichOpticalElement() << std::endl;*/
	//int b;std::cin >> b;
/*}*/
	/*if(icase == 2704 && iarg == BeforeTransport)
	{
		m_nDebugCounter++;
		std::cout << "-------------------" << m_nDebugCounter << "-----------------------" << std::endl;
		std::cout << "iarg:" << iarg << ",np-1:" << nStackAddress << ",ir:" << the_stack->ir[nStackAddress] - 2 << std::endl;
		std::cout << "r=(" << the_stack->x[nStackAddress] << "," << the_stack->y[nStackAddress]  << "," << the_stack->z[nStackAddress] <<")" << std::endl;
	}*/
	/*if(m_nDebugCounter >= 4505541)
	{
		std::cout << "-------------------" << m_nDebugCounter << "-----------------------" << std::endl;
		std::cout << "iarg:" << iarg << ",np-1:" << nStackAddress << ",ir:" << the_stack->ir[nStackAddress] - 2 << std::endl;
		std::cout << "r=(" << the_stack->x[nStackAddress] << "," << the_stack->y[nStackAddress]  << "," << the_stack->z[nStackAddress] <<")" << std::endl;
		std::cout << "p=(" << the_stack->u[nStackAddress] << "," << the_stack->v[nStackAddress]  << "," << the_stack->w[nStackAddress] <<")" << std::endl;
		std::cout << "nCharge: " << nCharge << std::endl;
		int aa;
		std::cin >> aa;
	}*/
	if(iarg == UserDiscard || iarg == PegsCut || iarg == EgsCut)
	{
		//reset all the parameters. i.e. check if propagated all paths produced by a grating
		m_pSplittingAlgorithm->PotentialParameterReset();
		e_fPhase[nStackAddress] = 0.0;
		e_fLogNorm[nStackAddress] = 0.0;
		e_bPrimary[nStackAddress] = false;
		the_stack->latch[nStackAddress] = 0;
	}
	else if(iarg == BeforeTransport && nCharge == 0/* && bIsPrimary == true*/)
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
			/*if(m_nDebugCounter >= 4505541)
			{
				std::cout << "tvstep: " << the_epcont->tvstep << std::endl;
			}*/
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
		//Splitting can only occur at medium interfaces
		if(nCurrentRegion != nPreviousRegion)
		{
			/*if(m_nDebugCounter >= 4505541)
			{
				std::cout << "before PotentialParticleSplitting" <<std::endl;
			}*/
			m_pSplittingAlgorithm->PotentialParticleSplitting();
			/*if(m_nDebugCounter >= 4505541)
			{
				std::cout << "after PotentialParticleSplitting" <<std::endl;
			}*/
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
					/*if(m_nDebugCounter >= 4505541)
					{
						std::cout << "before ApplySnellsLawHowfar" <<std::endl;
					}*/
					nReflectionTransmission = ApplySnellsLawHowfar(&fPosition, &m_fPositionBeforeTVStep, &fDirection, nCurrentRegion, nPreviousRegion,  nMediumIndexNew, nMediumIndexOld, fEnergy);
					/*if(m_nDebugCounter >= 4505541)
					{
						std::cout << "aftere ApplySnellsLawHowfar: " << nReflectionTransmission << std::endl;
					}*/
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
		if( iarg == BeforePhoto || iarg == BeforeCompton || iarg == BeforeRayleigh /*|| iarg == BeforePair*/)
		{
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

			//Check a few
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

void Example_Usercode::reportResults()
{
	egsInformation("\n\n last case = %d Etot = %g\n", (int)last_case,m_fEtot);
	std::cout << "Number of compton and photo effect scattering events of primary photons: " << m_nNumberOfPhotonScatteringEvents << std::endl;
	std::cout << "Number of Rayleigh scattering events of primary photons: " <<  m_nNumberOfRayleighEvents << std::endl;
	std::cout << "Pure ray tracing propagations: " << m_nNumberOfRayTracingPropagations << std::endl;
	std::cout << "Histories with too many reflections: " << m_nNumberOfParticlesWithToManyReflections <<std::endl;
	m_pSplittingAlgorithm->ReportSplittingSummary((int)last_case);
}


int main(int argc, char **argv)
{
	Example_Usercode simulation(argc, argv);

	simulation.run();

	simulation.reportResults();

	simulation.finish();

	return 0;
}
