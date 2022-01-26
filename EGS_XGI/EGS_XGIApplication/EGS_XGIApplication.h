#ifndef _EGS_XGIAPPLICATION_H_
#define _EGS_XGIAPPLICATION_H_

//General includes
#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>

//EGS includes
#include "egs_advanced_application.h"
#include "egs_interface2.h"
#include "egs_vector.h"
#include "egs_transformations.h"
#include "egs_interpolator.h"

#include "egs_functions.h"	/**/
#include "egs_base_source.h"	/**/
#include "egs_timer.h"		/**/
#include "egs_rndm.h"

#include "SplittingAlgorithm.h"
#include "RefractiveIndexCalculator.h"


class EGS_XGIApplication : public EGS_AdvancedApplication
{
public:
	EGS_XGIApplication(int argc, char **argv);

	virtual ~EGS_XGIApplication();

	virtual int ausgab(int iarg);

	virtual int run();

	virtual void startHistory(EGS_I64 this_case);

	virtual void endHistory();

	virtual void reportResults();

	virtual void finish();

	//Pure Ray Tracing (no interactions)
	int PerformRayTracing();


protected:

	void ActivateAusgabCals();

	EGS_I64	m_nNumberOfHistoriesToSimulate;
	EGS_I64 m_nReportProgress;
	EGS_Float m_fEtot;
	EGS_Float m_fNormalization;

	std::string m_sInputFileName;
	std::string m_sRefractiveIndexFileName;

	std::vector<EGS_Float> m_fMu_over_two_ForInitialEnergy;


	//A function to import refractive index deltas from the i_sRefFile file
	void ImportRefractiveIndex(string i_sRefFile);

	//i_pSurfaceNormal can be used to pass the surface normal for Snell's law
	//if i_pSurfaceNormal == 0 use howfar(.) to obtain the surface normal
	/*return value distinguishes 4 cases
		-1: Failed attempt
		0: refraction (transmission)
		1: reflection
		2: direction and surface normal are either orthogonal or parallel (up to an epsilon)
				i.e. no changes to the_stack
	*/
	int ApplySnellsLawHowfar(EGS_Vector* i_pPosition, EGS_Vector* i_pPositionOld, EGS_Vector* i_pMomentum, int i_nGlobalRegion, int i_nGlobalRegionOld, int i_nMediumIndex, int i_nMediumIndexOld, EGS_Float i_fPhotonEnergy, EGS_Vector* i_pSurfaceNormal = 0);

	unsigned int m_nNumberOfReflections;
	unsigned int m_nNumberOfParticlesWithToManyReflections;
	unsigned int m_nMaxNumberOfReflections;
	unsigned int m_nNumberOfReflectionsTotal;
	EGS_Vector m_fPositionBeforeTVStep;

	//Path splitting and scoring
	SplittingAlgorithm* m_pSplittingAlgorithm;

	int InitSplittingAlgorithm();
	//Try getting the refractive index and the attenuation coefficients, which are constant for one history in the current implementation

	int m_nNumberOfPhotonScatteringEvents;

	RefractiveIndexCalculator* m_pRefractiveIndexCalculator;
	int icase;
	EGS_Float m_fEnergyLastHistory;
	//int m_nDebugCounter;
private:
};


#endif
