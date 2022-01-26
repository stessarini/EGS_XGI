#ifndef _EXAMPLE_EXAMPLE_USERCODE_SCORE_ENERGY_H_
#define _EXAMPLE_EXAMPLE_USERCODE_SCORE_ENERGY_H_

#include <iostream>
#include <string.h>

#include "EGS_XGIApplication.h"
#include "egs_advanced_application.h"
#include "egs_interface2.h"

#include "egs_functions.h"	/**/
#include "egs_base_source.h"	/**/
#include "egs_timer.h"		/**/
#include "egs_rndm.h"
#include "egs_scoring.h"

class Example_Usercode_Score_Energy : public EGS_XGIApplication
{
private:

public:
	Example_Usercode_Score_Energy(int argc, char **argv);

	virtual ~Example_Usercode_Score_Energy();

	virtual int ausgab(int iarg);

	//virtual int run();

	//void ActivateAusgabCals();

	void startHistory(EGS_I64 this_case);

	void endHistory();

	void reportResults();

	//void finish();
protected:

	int m_nNumberOfRayTracingPropagations;
	int m_nNumberOfRayleighEvents;
//	int m_nStackSizeLastCall;
private:
	EGS_ScoringArray* m_pDepositedEnergy;
	int m_nNumberOfRegions;
	EGS_Float m_fScoredEnergyThisHistory;
	int m_nInteractionCheck;

};


#endif
