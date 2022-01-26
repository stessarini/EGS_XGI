#include <vector>
#include <iomanip>
#include <map>

#include "EGS_XGIApplication.h"
#include "DefaultSplittingAlgorithm.h"


//#include "gimc_stacks.h"
#include "xgi_global_variables.h"

#include "egs_functions.h"
#include "egs_interface2.h"
#include "egs_run_control.h"
#include "egs_base_source.h"
#include "egs_input.h"
#include "egs_interpolator.h"
#include "egs_rndm.h"
#include "egs_math.h"
#include "egs_ausgab_object.h"


EGS_XGIApplication::EGS_XGIApplication(int argc, char **argv)
	:EGS_AdvancedApplication(argc, argv),
	m_nNumberOfHistoriesToSimulate(1),
	m_nReportProgress(10),
	m_nNumberOfReflections(0),
	m_nNumberOfParticlesWithToManyReflections(0),
	m_nMaxNumberOfReflections(10),/*Just set it to a default value for now*/
	m_nNumberOfReflectionsTotal(0),
	m_fPositionBeforeTVStep(),
	m_pSplittingAlgorithm(NULL),
	m_fEtot(0.0),
	m_pRefractiveIndexCalculator(0)
{
	//Initialize some pointers
	/*i_ededx = NULL;
	i_pdedx = NULL;
	i_esig = NULL;
	i_psig = NULL;
	i_ebr1 = NULL;
	i_pbr1 = NULL;
	i_pbr2 = NULL;
	i_gmfp = NULL;
	i_gbr1 = NULL;
	i_gbr2 = NULL;
	i_cohe = NULL;
	i_photonuc = NULL;*/
	//SETUP MC
	m_nNumberOfPhotonScatteringEvents = 0;

	initSimulation();

	describeSimulation();

	ActivateAusgabCals();

    //Get parameters from argv
	for(int nArgumentIndex = 1; nArgumentIndex != argc; nArgumentIndex++)
	{
		if((strcmp("-i",argv[nArgumentIndex]) == 0))
		{
			m_sInputFileName = argv[nArgumentIndex + 1];

			std::cout << "The egs inputfile: " << m_sInputFileName << ".egsinp" << std::endl;
		}
		else if((strcmp("-r",argv[nArgumentIndex]) == 0))
		{
			m_sRefractiveIndexFileName = argv[nArgumentIndex + 1];
			std::cout << "The refractive index file: " << m_sRefractiveIndexFileName << ".refind" << std::endl;
		}
		else if((strcmp("-n",argv[nArgumentIndex]) == 0))
		{
			m_nNumberOfHistoriesToSimulate = stoll(argv[nArgumentIndex + 1]);
		}
		else if((strcmp("-nrep",argv[nArgumentIndex]) == 0))
		{
			m_nReportProgress = stoll(argv[nArgumentIndex + 1]);
		}
	}

	string sInputFile = m_sInputFileName;
	sInputFile += ".egsinp";

	//Import refractive index data
	string sRefractiveIndex = m_sRefractiveIndexFileName;
	sRefractiveIndex += ".refind";

	m_pRefractiveIndexCalculator = new RefractiveIndexCalculator();
	m_pRefractiveIndexCalculator->ImportRefractiveIndex(sRefractiveIndex, geometry);

	int errSpAlg = InitSplittingAlgorithm();//TODO avoid the start of a MC run without proper inputfile
	m_pSplittingAlgorithm->ReportSplittingAlgorithm();


	std::vector<EGS_Float> aux(geometry->nMedia(),0.0);
	m_fMu_over_two_ForInitialEnergy = aux;
}

EGS_XGIApplication::~EGS_XGIApplication()
{
	if (m_pRefractiveIndexCalculator != 0)
	{
		delete m_pRefractiveIndexCalculator;
	}
	if(m_pSplittingAlgorithm != 0)
	{
		delete m_pSplittingAlgorithm;
		m_pSplittingAlgorithm = 0;
	}
	delete [] i_ededx;
	delete [] i_pdedx;
	delete [] i_esig;
	delete [] i_psig;
	delete [] i_ebr1;
	delete [] i_pbr1;
	delete [] i_pbr2;
	delete [] i_gmfp;
	delete [] i_gbr1;
	delete [] i_gbr2;
	delete [] i_cohe;
	delete [] i_photonuc;

	//delete source;
	//delete geometry;
}

int EGS_XGIApplication::ausgab(int iarg)
{
	return 0;
}

int EGS_XGIApplication::run()
{
	egsInformation("\n\nSimulating %d particles...\n\n",m_nNumberOfHistoriesToSimulate);
  EGS_Timer timer;
  timer.start();
  EGS_I64 nperb = m_nNumberOfHistoriesToSimulate/m_nReportProgress;

	//reduce the number of log calculations
	m_fEnergyLastHistory = 0.0;
  if (nperb < 1)
  {
      nperb = 1;
  }
  EGS_Float aperb = ((EGS_Float)nperb)/((EGS_Float)m_nNumberOfHistoriesToSimulate);
  int q, latch;
  EGS_Float E, wt;
  EGS_Vector x,u;

  m_fNormalization = 0.0;

  int nNumberOfParticlesCreatedOutsideFOV = 0;
	unsigned int nNumberOfSimulatedHistories = 0;

  for (/*EGS_I64*/ icase=1; icase<=m_nNumberOfHistoriesToSimulate; icase++)
  {
    if (icase%nperb == 0)
    {
    	egsInformation("+Finished %7.2f percent of"
                                             " cases, cpu time = %9.3f\n",
                                             100*aperb*(icase/nperb),timer.time());
		}
		//Draw the next particle from the source.
		//rndm is a pointer to a random number generators.
    EGS_I64 this_case = source->getNextParticle(rndm,q,latch,E,wt,x,u);
    int ireg = geometry->isWhere(x);
    if (ireg < 0)
    {
      EGS_Float t = 1e30;
      ireg = geometry->howfar(ireg,x,u,t);
      nNumberOfParticlesCreatedOutsideFOV++;
      if (ireg >= 0)
      {
        x += u*t;
      }
    }
    if (ireg >= 0)
    {
			last_case = this_case;
     	if (q == 1)
      {
        m_fEtot += (E + 2*the_useful->rm)*wt;
        m_fNormalization += wt;
      }
      else
      {
        m_fEtot += E*wt;
        m_fNormalization += wt;
      }
			/*if(icase >= 2500)
			{
				std::cout << "icase: " << icase << std::endl;
			}
			if(icase == 2704)
			{
				std::cout << "case 2704" << std::endl;
				int a;
				std::cin >> a;
			}*/
			//m_nDebugCounter = 0;
			m_nNumberOfReflections = 0;
      the_stack->E[0] = (q) ? E + the_useful->rm : E;
      the_stack->x[0] = x.x;
      the_stack->y[0] = x.y;
      the_stack->z[0] = x.z;
      the_stack->u[0] = u.x;
      the_stack->v[0] = u.y;
      the_stack->w[0] = u.z;
      the_stack->dnear[0] = 0;
      the_stack->wt[0] = wt;
      the_stack->ir[0] = ireg+2;
      the_stack->iq[0] = q;
      the_stack->latch[0] = 0;//latch;
      the_stack->np = 1;
      e_fPhase[0] = 0.0;
      e_fLogNorm[0] = 0.0;
      e_bPrimary[0] = true;
			nNumberOfSimulatedHistories++;
			if(m_pSplittingAlgorithm->IsInInitialCondition() != 0)
			{
				std::cout << icase << ": " << m_pSplittingAlgorithm->WhichOpticalElement() << std::endl;
				std::cout << "Splitting algorithm not in initial state - Abort simulation" << std::endl;
				break;
			}
    	startHistory(this_case);
			m_fEnergyLastHistory = the_stack->E[0];
      egsShower();
      endHistory();
			m_nNumberOfReflectionsTotal += m_nNumberOfReflections;
    }
  }
  egsInformation("\n\nFinished simulation,CPU time was %g\n\n", timer.time());
  std::cout << "nNumberOfParticlesCreatedOutsideFOV = " << nNumberOfParticlesCreatedOutsideFOV << std::endl;
  std::cout << "m_nNumberOfPhotonScatteringEvents = " << m_nNumberOfPhotonScatteringEvents << std::endl;
	std::cout << "m_nNumberOfReflectionsTotal: " << m_nNumberOfReflectionsTotal << std::endl;
	std::cout << "m_nNumberOfParticlesWithToManyReflections: " << m_nNumberOfParticlesWithToManyReflections << std::endl;
  std::cout << "Average energy = " << m_fEtot / m_fNormalization << std::endl;
	m_pSplittingAlgorithm->EndSimulation(nNumberOfSimulatedHistories);
  return 0;
}

void EGS_XGIApplication::startHistory(EGS_I64 this_case)
{
	if(m_fEnergyLastHistory != the_stack->E[0])
	{
		m_pSplittingAlgorithm->StartRayTracing(the_stack->E[0]);
		for(unsigned int i = 0; i < geometry->nMedia(); i++)
		{
			m_fMu_over_two_ForInitialEnergy[i] = 1.0 / i_gmfp[i].interpolate(log(the_stack->E[0])) / 2.0;
		}
	}
}

void EGS_XGIApplication::endHistory()
{
	m_pSplittingAlgorithm->endHistory();
}

void EGS_XGIApplication::reportResults()
{
	egsInformation("\n\n last case = %d Etot = %g\n", (int)last_case,m_fEtot);
	m_pSplittingAlgorithm->ReportSplittingSummary((int)last_case);
}

void EGS_XGIApplication::finish()
{
	egsFinish();
}

void EGS_XGIApplication::ActivateAusgabCals()
{
	the_epcont->iausfl[BeforeTransport] = 1;     //!< before the step							0
	the_epcont->iausfl[EgsCut] = 1;              //!< energy below Ecut or Pcut					1
	the_epcont->iausfl[PegsCut] = 1;             //!< energy below AE or AP						2
	the_epcont->iausfl[UserDiscard ]= 1;         //!< user requested discard					3
	the_epcont->iausfl[ExtraEnergy] = 0;         /*!< initiated when part of the energy is not	4
														transfered to particles (e.g. binding energy)*/
	the_epcont->iausfl[AfterTransport] = 1;      //!< after the step							5
	the_epcont->iausfl[BeforeBrems] = 0;         //!< before a bremsstrahlung interaction		6
	the_epcont->iausfl[AfterBrems] = 0;          //!< after a bremsstrahlung interaction		7
	the_epcont->iausfl[BeforeMoller] = 0;        //!< before an inelastic collision (e-)		8
	the_epcont->iausfl[AfterMoller] = 0;         //!< after an inelastic collision (e-)			9
	the_epcont->iausfl[BeforeBhabha] = 0;       //!< before an inelastic collision (e+)			10
	the_epcont->iausfl[AfterBhabha] = 0;       //!< after an inelastic collision (e+)			11
	the_epcont->iausfl[BeforeAnnihFlight] = 0;  //!< before annihilation in flight				12
	the_epcont->iausfl[AfterAnnihFlight] = 0;   //!< after annihilation in flight				13
	the_epcont->iausfl[BeforeAnnihRest] = 0;    //!< before annihilation at rest				14 			28
	the_epcont->iausfl[AfterAnnihRest] = 0;     //!< after annihilation at rest					15			14
	the_epcont->iausfl[BeforePair] = 1;         //!< before pair production						16				15
	the_epcont->iausfl[AfterPair] = 1;          //!< after pair production						17				16
	the_epcont->iausfl[BeforeCompton] = 1;      //!< before a Compton scattering event			18	17
	the_epcont->iausfl[AfterCompton] = 1;       //!< after a Compton scattering event			19		18
	the_epcont->iausfl[BeforePhoto] = 1;        //!< before a photo-absorption event			20		19
	the_epcont->iausfl[AfterPhoto] = 1;         //!< after a photo-absorption event				21		20
	the_epcont->iausfl[EnteringUphi] = 0;       //!< the rotation routine was just entered		22	21
	the_epcont->iausfl[LeavingUphi] = 0;        //!< about to leave the rotation routine		23		22
		 //!< I consider the above 2 to be obsolete
	the_epcont->iausfl[BeforeRayleigh] = 1;     //!< before coherent scattering					24				23
	the_epcont->iausfl[AfterRayleigh] = 1;      //!< after coherent scattering					25				24
	the_epcont->iausfl[FluorescentEvent] = 0;   //!< a fluorescent transition just occured		26	25
	the_epcont->iausfl[CosterKronigEvent] = 0;  //!< a Coster-Kronig transition just occured	27	26
	the_epcont->iausfl[AugerEvent] = 0;         //!< an Auger transition just occured			28			27
	the_epcont->iausfl[BeforePhotoNuc] = 0;      //!< before a photonuclear event				29				29
	the_epcont->iausfl[AfterPhotoNuc] = 0;       //!< after a photonuclear event				30
	//the_epcont->iausfl[AfterSubPhoton] = 1;			//33 after sub-threshold photon energy deposition
	the_epcont->iausfl[UnknownCall] = 0;         //!< last element in the enumeration			31			35
}



/*This function performs ray-tracing using the top parameters on the_stack (e.g. the_stack->x[the_stack->np-1]). It transports the_stack entry untill it exits the geeometry.
IMPORTANT, this function changes the_stack, i.e. parameters on top of the_stack have to be saved in the meantime for later use*/
/*void EGS_XGIApplication::PerformRayTracing()
{
	//Get the initial parameters for the transport
	int nStackAddress = the_stack->np-1;
	//the_stack->latch[nStackAddress] = 2;
	int StacksizeBefore = nStackAddress;
	EGS_Float fEnergy = the_stack->E[nStackAddress];

	//Since energy is constant during ray-tracing
	std::vector<EGS_Float> fRefractiveIndex(geometry->nMedia());
	for(unsigned int nMediumCounter = 0; nMediumCounter != m_fMu_over_two_ForInitialEnergy.size(); nMediumCounter++)
	{
		fRefractiveIndex[nMediumCounter] = 1.0 - m_pRefractiveIndexCalculator->GetRefractiveIndexDelta(nMediumCounter, fEnergy);
	}

	//Get initial parameters from the_stack
	EGS_Vector fDirection;
	fDirection.x = the_stack->u[nStackAddress];
	fDirection.y = the_stack->v[nStackAddress];
	fDirection.z = the_stack->w[nStackAddress];
	fDirection.normalize();
	//in case u,v,w are copied during splitting:
	the_stack->u[nStackAddress] = fDirection.x;
	the_stack->v[nStackAddress] = fDirection.y;
	the_stack->w[nStackAddress] = fDirection.z;

	EGS_Vector fPositionOld;
	fPositionOld.x = the_stack->x[nStackAddress];
	fPositionOld.y = the_stack->y[nStackAddress];
	fPositionOld.z = the_stack->z[nStackAddress];

	EGS_Vector fPositionNew;

	//Overestimate the distance to the next interface
	EGS_Float fDistanceToNextInterface = 1e30;

	int nRegionIndexOld = the_stack->ir[nStackAddress] - 2;
	int nMediumOld = geometry->medium(nRegionIndexOld);
	int nMediumNew = -1;
	EGS_Vector fSurfaceNormal;
	int nRegionIndexNew= geometry->howfar(nRegionIndexOld, fPositionOld, fDirection, fDistanceToNextInterface, &nMediumNew, &fSurfaceNormal);
	//Don't do the +2 for region numbers, since don't use FORTRAN and don't do splitting based on region numbers.
	the_stack->ir[nStackAddress] = nRegionIndexNew;

	//Transport
	fPositionNew.x = fPositionOld.x + fDistanceToNextInterface * fDirection.x;
	fPositionNew.y = fPositionOld.y + fDistanceToNextInterface * fDirection.y;
	fPositionNew.z = fPositionOld.z + fDistanceToNextInterface * fDirection.z;

	the_stack->x[nStackAddress] = fPositionNew.x;
	the_stack->y[nStackAddress] = fPositionNew.y;
	the_stack->z[nStackAddress] = fPositionNew.z;

	//Track phase and norm
	if(nMediumOld >= 0)
	{
		e_fLogNorm[nStackAddress] -= fDistanceToNextInterface * m_fMu_over_two_ForInitialEnergy[nMediumOld];
		e_fPhase[nStackAddress] += fDistanceToNextInterface * fRefractiveIndex[nMediumOld];
	}
	else
	{
		//In case of vacuum vacuum
		e_fPhase[nStackAddress] += fDistanceToNextInterface;
	}

	//Split path
	m_pSplittingAlgorithm->PotentialParticleSplitting();

	int StackAddressAfterAppliedSplitting = the_stack->np -1;
	nRegionIndexNew = the_stack->ir[StackAddressAfterAppliedSplitting];
	fPositionNew.z = the_stack->z[StackAddressAfterAppliedSplitting];

	//if didn't split, new medium  and not leaving geometry, apply Snell's law
	if((StacksizeBefore == StackAddressAfterAppliedSplitting) && nMediumOld != nMediumNew && nRegionIndexNew != -1)
	{
		int nReflectionTransmission;
		nReflectionTransmission = ApplySnellsLawHowfar(&fPositionNew, &fPositionOld, &fDirection, nRegionIndexNew, nRegionIndexOld,  nMediumNew, nMediumOld, fEnergy, &fSurfaceNormal);
		if(abs(fDirection.length() - 1) > 1e-6 )
		{
			std::cout << "fDirection not normalized (" << fDirection.x <<", " << fDirection.y << "; " << fDirection.z << ")" << std::endl;
		}
		if(nReflectionTransmission == 0)
		{
			the_stack->u[nStackAddress] = fDirection.x;
			the_stack->v[nStackAddress] = fDirection.y;
			the_stack->w[nStackAddress] = fDirection.z;
		}
		else if (nReflectionTransmission == 1)
		{
			the_stack->u[nStackAddress] = fDirection.x;
			the_stack->v[nStackAddress] = fDirection.y;
			the_stack->w[nStackAddress] = fDirection.z;
			nRegionIndexNew = the_stack->ir[nStackAddress] - 2;
			nMediumNew = nMediumOld;
		}
		else if(nReflectionTransmission == -2)
		{
			std::cout << "Warning: Unknown call ApplySnellsLawHowfar" << std::endl;
		}
	}
	else if(nStackAddress < StackAddressAfterAppliedSplitting)
	{
		fDirection.x = the_stack->u[StackAddressAfterAppliedSplitting];
		fDirection.y = the_stack->v[StackAddressAfterAppliedSplitting];
		fDirection.z = the_stack->w[StackAddressAfterAppliedSplitting];
	}
	//the_stack->latch[nStackAddress] = 3;
	while( ((StacksizeBefore < the_stack->np-1) || nRegionIndexNew != -1) )
	{
		//Set the old values to new values
		fPositionOld = fPositionNew;
		fDistanceToNextInterface = 1e30;
		nStackAddress = the_stack->np-1;
		nRegionIndexOld = nRegionIndexNew;
		nMediumOld = nMediumNew;

		//Get new values
		nRegionIndexNew = geometry->howfar(nRegionIndexOld, fPositionOld, fDirection, fDistanceToNextInterface, &nMediumNew, &fSurfaceNormal);
		the_stack->ir[nStackAddress] = nRegionIndexNew;

		//Transport
		fPositionNew.x = fPositionOld.x + fDistanceToNextInterface * fDirection.x;
		fPositionNew.y = fPositionOld.y + fDistanceToNextInterface * fDirection.y;
		fPositionNew.z = fPositionOld.z + fDistanceToNextInterface * fDirection.z;

		the_stack->x[nStackAddress] = fPositionNew.x;
		the_stack->y[nStackAddress] = fPositionNew.y;
		the_stack->z[nStackAddress] = fPositionNew.z;

		//Track phase and norm
		if(nMediumOld >= 0)
		{
			e_fLogNorm[nStackAddress] -= fDistanceToNextInterface * m_fMu_over_two_ForInitialEnergy[nMediumOld];
			e_fPhase[nStackAddress] += fDistanceToNextInterface * fRefractiveIndex[nMediumOld];
		}
		else
		{
			e_fPhase[nStackAddress] += fDistanceToNextInterface;
		}

		m_pSplittingAlgorithm->PotentialParticleSplitting();
		StackAddressAfterAppliedSplitting = the_stack->np-1;
		//Update region number and position in case the splitting plane "absorbed" path
		nRegionIndexNew = the_stack->ir[StackAddressAfterAppliedSplitting];
		fPositionNew.z = the_stack->z[StackAddressAfterAppliedSplitting];
		if((nStackAddress == StackAddressAfterAppliedSplitting) && nMediumOld != nMediumNew && nRegionIndexNew != -1 )
		{
			int nReflectionTransmission;
			nReflectionTransmission =ApplySnellsLawHowfar(&fPositionNew, &fPositionOld, &fDirection, nRegionIndexNew, nRegionIndexOld,  nMediumNew, nMediumOld, fEnergy, &fSurfaceNormal);
			if(nReflectionTransmission == 0)
			{
				the_stack->u[nStackAddress] = fDirection.x;
				the_stack->v[nStackAddress] = fDirection.y;
				the_stack->w[nStackAddress] = fDirection.z;
			}
			else if (nReflectionTransmission == 1)
			{
				the_stack->u[nStackAddress] = fDirection.x;
				the_stack->v[nStackAddress] = fDirection.y;
				the_stack->w[nStackAddress] = fDirection.z;
				nRegionIndexNew = the_stack->ir[nStackAddress] - 2;
				nMediumNew = nMediumOld;
			}
			else if (nReflectionTransmission == -2)
			{
				std::cout << "Warming: Unknown call of ApplySnellsLawHowfar" << std::endl;
			}
		}
		else if(nStackAddress < StackAddressAfterAppliedSplitting)
		{
			fDirection.x = the_stack->u[StackAddressAfterAppliedSplitting];
			fDirection.y = the_stack->v[StackAddressAfterAppliedSplitting];
			fDirection.z = the_stack->w[StackAddressAfterAppliedSplitting];
		}
		if(nRegionIndexNew == -1 && (StacksizeBefore < StackAddressAfterAppliedSplitting))
		{
			the_stack->E[StackAddressAfterAppliedSplitting] = 0.0;

			the_stack->x[StackAddressAfterAppliedSplitting] = 0.0;
      the_stack->y[StackAddressAfterAppliedSplitting] = 0.0;
      the_stack->z[StackAddressAfterAppliedSplitting] = 0.0;

			the_stack->u[StackAddressAfterAppliedSplitting] = 0.0;
      the_stack->v[StackAddressAfterAppliedSplitting] = 0.0;
      the_stack->w[StackAddressAfterAppliedSplitting] = 0.0;

			the_stack->dnear[StackAddressAfterAppliedSplitting] = 0;
      the_stack->wt[StackAddressAfterAppliedSplitting] = 0.0;
      the_stack->ir[StackAddressAfterAppliedSplitting] = 0;
      the_stack->iq[StackAddressAfterAppliedSplitting] = 0;
      the_stack->latch[StackAddressAfterAppliedSplitting] = 0;

			e_fPhase[StackAddressAfterAppliedSplitting] = 0.0;
      e_fLogNorm[StackAddressAfterAppliedSplitting] = 0.0;
      e_bPrimary[StackAddressAfterAppliedSplitting] = false;

			nRegionIndexNew = the_stack->ir[StackAddressAfterAppliedSplitting - 1];

      fPositionNew.x = the_stack->x[StackAddressAfterAppliedSplitting -1];
      fPositionNew.y = the_stack->y[StackAddressAfterAppliedSplitting -1];
      fPositionNew.z = the_stack->z[StackAddressAfterAppliedSplitting -1];

			fDirection.x = the_stack->u[StackAddressAfterAppliedSplitting - 1];
			fDirection.y = the_stack->v[StackAddressAfterAppliedSplitting -1];
			fDirection.z = the_stack->w[StackAddressAfterAppliedSplitting - 1];

			nMediumNew = geometry->medium(nRegionIndexNew);

			m_pSplittingAlgorithm->PotentialParameterReset();
      the_stack->np -= 1;
		}
	}
	m_pSplittingAlgorithm->PotentialParameterReset();
	e_fPhase[StacksizeBefore] = 0.0;
	e_fLogNorm[StacksizeBefore] = 0.0;
}*/

int EGS_XGIApplication::PerformRayTracing()
{
	//int dummy;
	//Since energy is constant during ray-tracing
	//std::cout << "RT starts: " << the_stack->np-1 << std::endl;
	int nStackAddressBeforeRayTracing = the_stack->np-1;
	int nRegionIndexOld = the_stack->ir[nStackAddressBeforeRayTracing] - 2;
	/*if(m_nDebugCounter >= 4505541)
	{
		std::cout << "-------------------" << "RT start" << "-----------------------" << std::endl;
		std::cout << ",np-1:" << nStackAddressBeforeRayTracing << ",ir:" << the_stack->ir[nStackAddressBeforeRayTracing] - 2 << std::endl;
		std::cout << "r=(" << the_stack->x[nStackAddressBeforeRayTracing] << "," << the_stack->y[nStackAddressBeforeRayTracing]  << "," << the_stack->z[nStackAddressBeforeRayTracing] <<")" << std::endl;
	}*/
	if(nRegionIndexOld < 0)
	{
		std::cout << "RT initial position outside..." << std::endl;
		//UserDiscard
		m_pSplittingAlgorithm->PotentialParameterReset();
		//delete the path
		the_stack->E[nStackAddressBeforeRayTracing] = 0.0;

		the_stack->x[nStackAddressBeforeRayTracing] = 0.0;
		the_stack->y[nStackAddressBeforeRayTracing] = 0.0;
		the_stack->z[nStackAddressBeforeRayTracing] = 0.0;

		the_stack->u[nStackAddressBeforeRayTracing] = 0.0;
		the_stack->v[nStackAddressBeforeRayTracing] = 0.0;
		the_stack->w[nStackAddressBeforeRayTracing] = 0.0;

		the_stack->dnear[nStackAddressBeforeRayTracing] = 0;
		the_stack->wt[nStackAddressBeforeRayTracing] = 0.0;
		the_stack->ir[nStackAddressBeforeRayTracing] = 1;
		the_stack->iq[nStackAddressBeforeRayTracing] = 0;
		the_stack->latch[nStackAddressBeforeRayTracing] = 0;

		e_fPhase[nStackAddressBeforeRayTracing] = 0.0;
		e_fLogNorm[nStackAddressBeforeRayTracing] = 0.0;
		e_bPrimary[nStackAddressBeforeRayTracing] = false;
		the_stack->np = the_stack->np-1;
		return 0;
	}
  EGS_Float fEnergy = the_stack->E[nStackAddressBeforeRayTracing];
	std::vector<EGS_Float> fRefractiveIndex(geometry->nMedia());
	for(unsigned int nMediumCounter = 0; nMediumCounter != m_fMu_over_two_ForInitialEnergy.size(); nMediumCounter++)
	{
		fRefractiveIndex[nMediumCounter] = 1.0 - m_pRefractiveIndexCalculator->GetRefractiveIndexDelta(nMediumCounter, fEnergy);
	}

	//BeforeTransport:

  //Get parameters from the_stack
    /*
		position
    Direction
    region index
    */
	EGS_Vector fDirection;
	fDirection.x = the_stack->u[nStackAddressBeforeRayTracing];
	fDirection.y = the_stack->v[nStackAddressBeforeRayTracing];
	fDirection.z = the_stack->w[nStackAddressBeforeRayTracing];
	fDirection.normalize();
	//in case u,v,w are copied during splitting:
	the_stack->u[nStackAddressBeforeRayTracing] = fDirection.x;
	the_stack->v[nStackAddressBeforeRayTracing] = fDirection.y;
	the_stack->w[nStackAddressBeforeRayTracing] = fDirection.z;

	EGS_Vector fPositionOld;
	fPositionOld.x = the_stack->x[nStackAddressBeforeRayTracing];
	fPositionOld.y = the_stack->y[nStackAddressBeforeRayTracing];
	fPositionOld.z = the_stack->z[nStackAddressBeforeRayTracing];

  int nMediumOld = geometry->medium(nRegionIndexOld);
	int nMediumNew = -1;
	EGS_Vector fSurfaceNormal;
	//std::cout << ".............................................." << std::endl;
	//std::cout << "RT before x = " << fPositionOld.x << ", " << fPositionOld.y << ", " << fPositionOld.z << std::endl;
	//std::cout << "medium name: " << geometry->getMediumName(nMediumOld) << std::endl;
	//std::cout << "np-1 = " << the_stack->np-1 << std::endl;
	/*if(m_nDebugCounter >= 4505541)
	{
		std::cout << ",np-1:" << nStackAddressBeforeRayTracing << ",ir:" << the_stack->ir[nStackAddressBeforeRayTracing] - 2 << std::endl;
		std::cout << "r=(" << fPositionOld.x << "," << fPositionOld.y  << "," << fPositionOld.z <<")" << std::endl;
		std::cout << "nMediumOld: " << nMediumOld << std::endl;
	}*/

  EGS_Float fDistanceToNextInterface = 1e30;
  int nRegionIndexNew= geometry->howfar(nRegionIndexOld, fPositionOld, fDirection, fDistanceToNextInterface, &nMediumNew, &fSurfaceNormal);
	//track phase and norm
	if(nMediumOld >= 0)
	{
		e_fLogNorm[nStackAddressBeforeRayTracing] -= fDistanceToNextInterface * m_fMu_over_two_ForInitialEnergy[nMediumOld];
		e_fPhase[nStackAddressBeforeRayTracing] += fDistanceToNextInterface * fRefractiveIndex[nMediumOld];
	}
	else
	{
		//In case of vacuum
		e_fPhase[nStackAddressBeforeRayTracing] += fDistanceToNextInterface;
	}

	//Do the RT step
  EGS_Vector fPositionNew;
  the_stack->ir[nStackAddressBeforeRayTracing] = nRegionIndexNew + 2;
  fPositionNew.x = fPositionOld.x + fDistanceToNextInterface * fDirection.x;
	fPositionNew.y = fPositionOld.y + fDistanceToNextInterface * fDirection.y;
	fPositionNew.z = fPositionOld.z + fDistanceToNextInterface * fDirection.z;
  the_stack->x[nStackAddressBeforeRayTracing] = fPositionNew.x;
  the_stack->y[nStackAddressBeforeRayTracing] = fPositionNew.y;
  the_stack->z[nStackAddressBeforeRayTracing] = fPositionNew.z;

	/*if(m_nDebugCounter >= 4505541)
	{
		std::cout << ",np-1:" << nStackAddressBeforeRayTracing << ",ir:" << the_stack->ir[nStackAddressBeforeRayTracing] - 2 << std::endl;
		std::cout << "r_new=(" << fPositionNew.x << "," << fPositionNew.y  << "," << fPositionNew.z <<")" << std::endl;
		std::cout << "nMediumNew: " << nMediumNew << std::endl;
	}*/
	/*std::cout << "RT after x = " << fPositionNew.x << ", " << fPositionNew.y << ", " << fPositionNew.z << std::endl;
	if(nMediumNew >= 0)
	{
		std::cout << "medium name: " << geometry->getMediumName(nMediumNew) << std::endl;
	}
	std::cout << "np-1 = " << the_stack->np-1 << std::endl;*/
	//AfterTransport:

  //Split path
	m_pSplittingAlgorithm->PotentialParticleSplitting();
  int StackAddressAfterAppliedSplitting = the_stack->np -1;
	/*if(m_nDebugCounter >= 4505541)
	{
		std::cout << "after splitting: np-1:" << StackAddressAfterAppliedSplitting << ",ir:" << the_stack->ir[StackAddressAfterAppliedSplitting] - 2 << std::endl;
	}*/
  if(nStackAddressBeforeRayTracing == StackAddressAfterAppliedSplitting)
  {
  	nRegionIndexNew = the_stack->ir[nStackAddressBeforeRayTracing] - 2;
  	//fPositionNew.z = the_stack->z[nStackAddressBeforeRayTracing];
    if(nRegionIndexNew < 0)
    {
			//UserDiscard
			m_pSplittingAlgorithm->PotentialParameterReset();
			//delete the path
			the_stack->E[nStackAddressBeforeRayTracing] = 0.0;

			the_stack->x[nStackAddressBeforeRayTracing] = 0.0;
			the_stack->y[nStackAddressBeforeRayTracing] = 0.0;
			the_stack->z[nStackAddressBeforeRayTracing] = 0.0;

			the_stack->u[nStackAddressBeforeRayTracing] = 0.0;
			the_stack->v[nStackAddressBeforeRayTracing] = 0.0;
			the_stack->w[nStackAddressBeforeRayTracing] = 0.0;

			the_stack->dnear[nStackAddressBeforeRayTracing] = 0;
			the_stack->wt[nStackAddressBeforeRayTracing] = 0.0;
			the_stack->ir[nStackAddressBeforeRayTracing] = 1;
			the_stack->iq[nStackAddressBeforeRayTracing] = 0;
			the_stack->latch[nStackAddressBeforeRayTracing] = 0;

			e_fPhase[nStackAddressBeforeRayTracing] = 0.0;
			e_fLogNorm[nStackAddressBeforeRayTracing] = 0.0;
			e_bPrimary[nStackAddressBeforeRayTracing] = false;
			the_stack->np = the_stack->np-1;
      return 0;
    }
    else if(nMediumOld != nMediumNew)
    {
			/*if(m_nDebugCounter >= 4505541)
			{
				std::cout << "Before SnellsLaw" << std::endl;
			}*/
      int nReflectionTransmission = ApplySnellsLawHowfar(&fPositionNew, &fPositionOld, &fDirection, nRegionIndexNew, nRegionIndexOld,  nMediumNew, nMediumOld, fEnergy, &fSurfaceNormal);
			/*if(m_nDebugCounter >= 4505541)
			{
				std::cout << "after SnellsLaw" << std::endl;
			}*/
			if(nReflectionTransmission == 0)
  		{
  			the_stack->u[nStackAddressBeforeRayTracing] = fDirection.x;
  			the_stack->v[nStackAddressBeforeRayTracing] = fDirection.y;
  			the_stack->w[nStackAddressBeforeRayTracing] = fDirection.z;
  		}
  		else if (nReflectionTransmission == 1)
  		{
  			the_stack->u[nStackAddressBeforeRayTracing] = fDirection.x;
  			the_stack->v[nStackAddressBeforeRayTracing] = fDirection.y;
  			the_stack->w[nStackAddressBeforeRayTracing] = fDirection.z;
  			nRegionIndexNew = the_stack->ir[nStackAddressBeforeRayTracing] - 2;
  			nMediumNew = nMediumOld;
  		}
  		else if(nReflectionTransmission == -2)
  		{
  			std::cout << "Warning: Unknown call ApplySnellsLawHowfar" << std::endl;
  		}
    }
  }
	else //if(nStackAddress < StackAddressAfterAppliedSplitting)
  {
    fDirection.x = the_stack->u[StackAddressAfterAppliedSplitting];
		fDirection.y = the_stack->v[StackAddressAfterAppliedSplitting];
		fDirection.z = the_stack->w[StackAddressAfterAppliedSplitting];
  }
	//std::cin >> dummy;
  int nStackAddress = StackAddressAfterAppliedSplitting;
	/*if(m_nDebugCounter >= 4505541)
	{
		std::cout << "Before while loop RT" << std::endl;
	}*/
  while(nStackAddressBeforeRayTracing <= nStackAddress)
  {
    //BeforeTransport
    //set old values
    fPositionOld = fPositionNew;
    nRegionIndexOld = nRegionIndexNew;
    nMediumOld = nMediumNew;
		/*if(m_nDebugCounter >= 4505541)
		{
			std::cout << "++++++++++++++++++++++++++++++++" << std::endl;
			std::cout << ",np-1:" << nStackAddress << ",ir:" << the_stack->ir[nStackAddress] - 2 << std::endl;
			std::cout << "r=(" << fPositionOld.x << "," << fPositionOld.y  << "," << fPositionOld.z <<")" << std::endl;
			std::cout << "nMediumOld: " << nMediumOld << std::endl;
			int aa;
			std::cin >> aa;
		}*/
		/*std::cout << ".............................................." << std::endl;
		std::cout << "RT before x = " << fPositionOld.x << ", " << fPositionOld.y << ", " << fPositionOld.z << std::endl;
		std::cout << "medium name: " << geometry->getMediumName(nMediumOld) << std::endl;
		std::cout << "np-1 = " << the_stack->np-1 << std::endl;*/
    fDistanceToNextInterface = 1e30;
    nRegionIndexNew = geometry->howfar(nRegionIndexOld, fPositionOld, fDirection, fDistanceToNextInterface, &nMediumNew, &fSurfaceNormal);

		//Track phase and norm
		if(nMediumOld >= 0)
		{
			e_fLogNorm[nStackAddress] -= fDistanceToNextInterface * m_fMu_over_two_ForInitialEnergy[nMediumOld];
			e_fPhase[nStackAddress] += fDistanceToNextInterface * fRefractiveIndex[nMediumOld];
		}
		else
		{
			e_fPhase[nStackAddress] += fDistanceToNextInterface;
		}

    //Transport step
    the_stack->ir[nStackAddress] = nRegionIndexNew + 2;
		fPositionNew.x = fPositionOld.x + fDistanceToNextInterface * fDirection.x;
		fPositionNew.y = fPositionOld.y + fDistanceToNextInterface * fDirection.y;
		fPositionNew.z = fPositionOld.z + fDistanceToNextInterface * fDirection.z;
    the_stack->x[nStackAddress] = fPositionNew.x;
		the_stack->y[nStackAddress] = fPositionNew.y;
		the_stack->z[nStackAddress] = fPositionNew.z;
		/*if(m_nDebugCounter >= 4505541)
		{
			std::cout << ",np-1:" << nStackAddress << ",ir:" << the_stack->ir[nStackAddress] - 2 << std::endl;
			std::cout << "r_new=(" << fPositionNew.x << "," << fPositionNew.y  << "," << fPositionNew.z <<")" << std::endl;
			std::cout << "nMediumNew: " << nMediumNew << std::endl;
		}*/
		/*std::cout << "RT after x = " << fPositionNew.x << ", " << fPositionNew.y << ", " << fPositionNew.z << std::endl;
		std::cout << "medium name: " << geometry->getMediumName(nMediumNew) << std::endl;
		std::cout << "np-1 = " << the_stack->np-1 << std::endl;*/

    //AfterTransport
    m_pSplittingAlgorithm->PotentialParticleSplitting();
		StackAddressAfterAppliedSplitting = the_stack->np-1;
		/*if(m_nDebugCounter >= 4505541)
		{
			std::cout << "after splitting: np-1:" << StackAddressAfterAppliedSplitting << ",ir:" << the_stack->ir[StackAddressAfterAppliedSplitting] - 2 << std::endl;
		}*/
    if(nStackAddress == StackAddressAfterAppliedSplitting)
    {
    	nRegionIndexNew = the_stack->ir[nStackAddress] - 2;
    	//fPositionNew.z = the_stack->z[nStackAddress];
      if(nRegionIndexNew < 0)
      {
				//UserDiscard
				//std::cout << "RT Discard in while loop" << std::endl;
				m_pSplittingAlgorithm->PotentialParameterReset();
        //delete the path
        the_stack->E[StackAddressAfterAppliedSplitting] = 0.0;

  			the_stack->x[StackAddressAfterAppliedSplitting] = 0.0;
        the_stack->y[StackAddressAfterAppliedSplitting] = 0.0;
        the_stack->z[StackAddressAfterAppliedSplitting] = 0.0;

  			the_stack->u[StackAddressAfterAppliedSplitting] = 0.0;
        the_stack->v[StackAddressAfterAppliedSplitting] = 0.0;
        the_stack->w[StackAddressAfterAppliedSplitting] = 0.0;

  			the_stack->dnear[StackAddressAfterAppliedSplitting] = 0;
        the_stack->wt[StackAddressAfterAppliedSplitting] = 0.0;
        the_stack->ir[StackAddressAfterAppliedSplitting] = 1;
        the_stack->iq[StackAddressAfterAppliedSplitting] = 0;
        the_stack->latch[StackAddressAfterAppliedSplitting] = 0;

  			e_fPhase[StackAddressAfterAppliedSplitting] = 0.0;
        e_fLogNorm[StackAddressAfterAppliedSplitting] = 0.0;
        e_bPrimary[StackAddressAfterAppliedSplitting] = false;

				nStackAddress--;
        the_stack->np = the_stack->np - 1;

				//Get parameters of the next stack entry
  			nRegionIndexNew = the_stack->ir[nStackAddress] - 2;

        fPositionNew.x = the_stack->x[nStackAddress];
        fPositionNew.y = the_stack->y[nStackAddress];
        fPositionNew.z = the_stack->z[nStackAddress];

  			fDirection.x = the_stack->u[nStackAddress];
  			fDirection.y = the_stack->v[nStackAddress];
  			fDirection.z = the_stack->w[nStackAddress];

  			nMediumNew = geometry->medium(nRegionIndexNew);
				//std::cout << "np-1 = " << the_stack->np-1 << std::endl;

      }
      else if(nMediumOld != nMediumNew)
      {
				/*if(m_nDebugCounter >= 4505541)
				{
					std::cout << "before SnellsLaw" <<std::endl;
				}*/
        int nReflectionTransmission = ApplySnellsLawHowfar(&fPositionNew, &fPositionOld, &fDirection, nRegionIndexNew, nRegionIndexOld,  nMediumNew, nMediumOld, fEnergy, &fSurfaceNormal);
				/*if(m_nDebugCounter >= 4505541)
				{
					std::cout << "after SnellsLaw " << nReflectionTransmission <<std::endl;
				}*/
        if(nReflectionTransmission == 0)
    		{
    			the_stack->u[nStackAddress] = fDirection.x;
    			the_stack->v[nStackAddress] = fDirection.y;
    			the_stack->w[nStackAddress] = fDirection.z;
    		}
    		else if (nReflectionTransmission == 1)
    		{
    			the_stack->u[nStackAddress] = fDirection.x;
    			the_stack->v[nStackAddress] = fDirection.y;
    			the_stack->w[nStackAddress] = fDirection.z;
    			nRegionIndexNew = the_stack->ir[nStackAddress] - 2;
    			nMediumNew = nMediumOld;
    		}
    		else if(nReflectionTransmission == -2)
    		{
    			std::cout << "Warning: Unknown call ApplySnellsLawHowfar" << std::endl;
    		}
      }
    }
  	else
    {
      fDirection.x = the_stack->u[StackAddressAfterAppliedSplitting];
  		fDirection.y = the_stack->v[StackAddressAfterAppliedSplitting];
  		fDirection.z = the_stack->w[StackAddressAfterAppliedSplitting];
      nStackAddress = StackAddressAfterAppliedSplitting;
    }
  }
	/*std::cout << "RT ends: " << the_stack->np-1 << std::endl;
	std::cin >> dummy;*/
	/*if(m_nDebugCounter >= 4505541)
	{
		std::cout << "RT end" <<std::endl;
	}*/
  return 0;
}


int EGS_XGIApplication::ApplySnellsLawHowfar(EGS_Vector* i_pPosition, EGS_Vector* i_pPositionOld, EGS_Vector* i_pMomentum, int i_nGlobalRegion, int i_nGlobalRegionOld, int i_nMediumIndex, int i_nMediumIndexOld, EGS_Float i_fPhotonEnergy, EGS_Vector* i_pSurfaceNormal)
{
	EGS_Vector SurfaceNormal;
	int region;
	EGS_Float fDistanceToNextInterface = 1e30;
	int nMediumNew;
	int o_nReport = -2;
	//According to pirs701 p. 114 the momentum direction is not necessarily normalized.
	i_pMomentum->normalize();
	if(i_pSurfaceNormal != 0)
	{
		SurfaceNormal = *i_pSurfaceNormal;
	}
	else
	{
		region = geometry->howfar(i_nGlobalRegionOld, *i_pPositionOld, *i_pMomentum, fDistanceToNextInterface, &nMediumNew, &SurfaceNormal);
	}
	if(abs(SurfaceNormal.length() -1.0) >= 1e-6)
	{
		std::cout << "#####################################" << std::endl;
		cout << "Refraction::SnellsLaw::ERROR: The surface normal is not normalized" << endl;
		cout << "Global Regions: " << i_nGlobalRegion << ", " << i_nGlobalRegionOld << endl;
		cout << "position: (" << i_pPosition->x << ", " << i_pPosition->y << ", " << i_pPosition->z << ")" << endl;
		cout << "position_old: (" << i_pPositionOld->x << ", " << i_pPositionOld->y << ", " << i_pPositionOld->z << ")" << endl;
		cout << "direction: (" << i_pMomentum->x << ", " << i_pMomentum->y << ", " << i_pMomentum->z << ")" << endl;
		cout << "surface normal from howfar :" << SurfaceNormal.x << ", " << SurfaceNormal.y << ", " << SurfaceNormal.z << std::endl;
		cout << e_bPrimary[the_stack->np -1] << std::endl;
		cout << "fDistanceToNextInterface = " << fDistanceToNextInterface << std::endl;
		std::cout << "howfar returns: " << region << std::endl;

		std::cout << "howfar(" << i_nGlobalRegionOld << ", " << "(" << i_pPositionOld->x << ", " << i_pPositionOld->y << ", " << i_pPositionOld->z << "), (" << i_pMomentum->x << ", " << i_pMomentum->y << ", " << i_pMomentum->z << "), " << 1e30 << "...)" << std::endl;
		std::cout << "primary: ";
		if(e_bPrimary[the_stack->np-1] == true)
		{
			std::cout << "true" << std::endl;
		}
		else
		{
			std::cout << "false" << std::endl;
		}
		std::cout << "latch = " << the_stack->latch[the_stack->np-1] << std::endl;
		std::cout << "Reflections: " << m_nNumberOfReflections << std::endl;
		std::cout << "old medium name :" <<  geometry->getMediumName(i_nGlobalRegionOld) << std::endl;
		std::cout << "new medium name :" <<  geometry->getMediumName(i_nGlobalRegion) << std::endl;
		std::cout << "geom name: " << geometry->getName() << std::endl;
		o_nReport = -1;
	}
	else
	{
		//get the cosine of the angle between the SurfaceNormal and the momentum direction:
		EGS_Float fCosAlpha = (SurfaceNormal * (*i_pMomentum));
		EGS_Float fAbsCosAlpha = abs(fCosAlpha);
		//Only apply Snell's law if SurfaceNormal and momentum direction are not orthogonal or tangential
		if(fAbsCosAlpha > 1e-6 && (1.0 - fAbsCosAlpha) > 1e-6)
		{
			//Ensure that the surface normal vector points into the same region as the momentum direction
			if(fCosAlpha < 0.0)
			{
				SurfaceNormal  = SurfaceNormal * (-1);
				fCosAlpha = -fCosAlpha;
			}

			//Get the rotation axis for the Snells Law (cross product)
			EGS_Vector RotationAxis = (*i_pMomentum) % SurfaceNormal;
			//Get abs(sin(alpha)) (=fSinAlpha)
			EGS_Float fSinAlpha = RotationAxis.length();
			//Normalize the rotation axis
			RotationAxis.normalize();

			EGS_Float fRefractiveIndex = 1.0 - m_pRefractiveIndexCalculator->GetRefractiveIndexDelta(i_nMediumIndex, i_fPhotonEnergy);
			EGS_Float fRefractiveIndexOld = 1.0 - m_pRefractiveIndexCalculator->GetRefractiveIndexDelta(i_nMediumIndexOld, i_fPhotonEnergy);
			//Angles according to Snell's law
			EGS_Float fSinBeta = (fRefractiveIndexOld/fRefractiveIndex) * fSinAlpha;
			//the rotation angle for the rotation matrix
			EGS_Float fPhi = 0.0;

			//Distinguish between transmission and reflection
			if(abs(fSinBeta)>1)//total reflection
			{
				if(m_nNumberOfReflections < m_nMaxNumberOfReflections)
				{
					fPhi = 2 * asin(fSinAlpha) + ec_fPi;
					int nNumberOfParticlesOnStack = the_stack->np -1;
					the_stack->ir[nNumberOfParticlesOnStack] = i_nGlobalRegionOld + 2;
					the_epcont->irnew = i_nGlobalRegionOld + 2;
					m_nNumberOfReflections++;
					o_nReport = 1;
					e_fPhase[nNumberOfParticlesOnStack] += ec_fPi;
					if (m_nNumberOfReflections == m_nMaxNumberOfReflections)
					{
						m_nNumberOfParticlesWithToManyReflections++;
					}
				}
				else
				{
					//o_nReport = 3;
					return 3;
				}
			}
			else
			{
				fPhi = asin(fSinAlpha) - asin(fSinBeta);
				o_nReport = 0;
			}
			//Get the rotation matrix
			//define some auxillary variables:
			EGS_Float C = cos(fPhi);
			EGS_Float S = sin(fPhi);
			EGS_Float t = 1-C;
			//The rotation matrix for Snells law
			EGS_RotationMatrix SnellsLaw(	t * RotationAxis.x * RotationAxis.x + C, t*RotationAxis.x * RotationAxis.y - S* RotationAxis.z, t*RotationAxis.x * RotationAxis.z + S * RotationAxis.y,
							t * RotationAxis.x * RotationAxis.y + S * RotationAxis.z, t* RotationAxis.y * RotationAxis.y + C, t * RotationAxis.y * RotationAxis.z - S * RotationAxis.x,
							t * RotationAxis.x * RotationAxis.z - S * RotationAxis.y, t * RotationAxis.y * RotationAxis.z + S * RotationAxis.x, t * RotationAxis.z * RotationAxis.z + C);

			//Apply the rotation matrix to the momentum direction vector
			*i_pMomentum = SnellsLaw * (*i_pMomentum);
			i_pMomentum->normalize();
		}
		else
		{
			o_nReport = 2;
		}
	}
	return o_nReport;
}

int EGS_XGIApplication::InitSplittingAlgorithm()
{
	EGS_Input* pSplittingAlgDefInput = NULL;
	bool bdelete_input = false;
	//input->print(0,std::cout);
	if(!input->isA("splitting algorithm definition"))
	{
		pSplittingAlgDefInput = input->takeInputItem("splitting algorithm definition");
		bdelete_input = true;
	}
	if(!pSplittingAlgDefInput)
	{
		std::cout << "EGS_XGIApplication::InitSplittingAlgorithm: No splitting algorithm definition" << std::endl;
		return 0;
	}
	EGS_Input* pSplittingAlgInput = NULL;
	if(!pSplittingAlgDefInput->isA("splitting algorithm"))
	{
		pSplittingAlgInput = pSplittingAlgDefInput->takeInputItem("splitting algorithm");
	}
	if(pSplittingAlgDefInput == NULL)
	{
		std::cout << "EGS_XGIApplication::InitSplittingAlgorithm: No splitting algorithm section" << std::endl;
		return 0;
	}
	else
	{
		std::vector<std::string> sType;
		int errSplittingAlgorithmType = pSplittingAlgInput->getInput("type",sType);
		if(!errSplittingAlgorithmType)
		{
			if(sType[0] == "DefaultSplittingAlgorithm")
			{
				m_pSplittingAlgorithm = new DefaultSplittingAlgorithm();
			}
			else
			{
				std::cout << "EGS_XGIApplication::InitSplittingAlgorithm: Unknown splitting algotrithm type" << std::endl;
				return 0;
			}

			if(m_pSplittingAlgorithm != NULL)
			{
				int errSpAlg = m_pSplittingAlgorithm->InitSplittingAlgorithm(pSplittingAlgDefInput, &m_sInputFileName, rndm, geometry, i_gmfp, m_pRefractiveIndexCalculator);
				if(errSpAlg == 0)
				{
					std::cout << "EGS_XGIApplication::InitSplittingAlgorithm: Error during 'InitSplittingAlgorithm(.)'" << std::endl;
					return 0;
				}
				int errDerSpAlg = m_pSplittingAlgorithm->InitDerivedSplittingAlgorithm(pSplittingAlgInput);
				if(errDerSpAlg == 0)
				{
					std::cout << "EGS_XGIApplication::InitSplittingAlgorithm: Error during 'InitDerivedSplittingAlgorithm(.)'" << std::endl;
					return 0;
				}
			}
			delete pSplittingAlgInput;
			delete pSplittingAlgDefInput;
		}
		return 1;
	}
}
