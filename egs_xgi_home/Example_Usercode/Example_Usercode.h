#ifndef _EXAMPLE_Example_Usercode_H_
#define _EXAMPLE_Example_Usercode_H_

#include <iostream>
#include <string.h>

#include "EGS_XGIApplication.h"
#include "egs_advanced_application.h"
#include "egs_interface2.h"

#include "egs_functions.h"	/**/
#include "egs_base_source.h"	/**/
#include "egs_timer.h"		/**/
#include "egs_rndm.h"


class Example_Usercode : public EGS_XGIApplication
{
private:

public:
	Example_Usercode(int argc, char **argv);

	virtual ~Example_Usercode();

	virtual int ausgab(int iarg);

	//virtual int run();

	//void ActivateAusgabCals();

	//void startHistory(EGS_I64 this_case);

	//void endHistory();

	void reportResults();

	//void finish();
protected:

	int m_nNumberOfRayTracingPropagations;
	int m_nNumberOfRayleighEvents;
//	int m_nStackSizeLastCall;
private:

};


#endif
