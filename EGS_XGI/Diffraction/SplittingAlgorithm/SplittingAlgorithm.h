#ifndef _SPLITTINGALGORITHM_H_
#define _SPLITTINGALGORITHM_H_
#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>

#include "OpticalElement.h"
#include "RefractiveIndexCalculator.h"

#include "egs_rndm.h"
#include "egs_interface2.h"
#include "egs_input.h"


class SplittingAlgorithm
{
public:

SplittingAlgorithm();

virtual ~SplittingAlgorithm();

SplittingAlgorithm(std::string i_sInputFile);

//Report relevant variables of this and OpticalElements
virtual void ReportSplittingAlgorithm();

//Report of OpticalElements after run
virtual void ReportSplittingSummary(int ncase);

//Do a particle splitting if necessary
virtual	void PotentialParticleSplitting();

//Reset Parameters, e.g. pointer to OpticalElement
virtual	void PotentialParameterReset();

int InitSplittingAlgorithm(EGS_Input* i_Input, std::string* i_sInputFileName, EGS_RandomGenerator* i_RandomNumberGenerator, EGS_BaseGeometry* i_pGeometry, EGS_Interpolator* i_gmfp, RefractiveIndexCalculator* i_pRefractiveIndexCalculator);

virtual int InitDerivedSplittingAlgorithm(EGS_Input* i_Input);

OpticalElement* GetPointerTo(std::string i_sName);

virtual void StartRayTracing(EGS_Float i_fInitialEnergy);

virtual void endHistory();

virtual void EndSimulation(int i_nNumberOfHistories);

virtual std::string WhichOpticalElement();

//return 0 if the splitting algorithm is in the correct state to start a new history
//otherwise erturn non-zero int
virtual int IsInInitialCondition();

protected:
vector<OpticalElement*> m_pSplittingObjects;
vector<OpticalElement*> m_pAllSplittingObjectsInInputfile; //All objects defined in the input file

};


#endif
