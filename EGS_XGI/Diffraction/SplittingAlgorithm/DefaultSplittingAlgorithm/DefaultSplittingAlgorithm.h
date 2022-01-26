#ifndef _DEFAULTSPLITTINGALGORITHM_H_
#define _DEFAULTSPLITTINGALGORITHM_H_
#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>

#include "SplittingAlgorithm.h"
//#include "SplittingObject.h"
#include "egs_input.h"

#include "egs_rndm.h"



class DefaultSplittingAlgorithm : public SplittingAlgorithm
{
public:
	DefaultSplittingAlgorithm();

	~DefaultSplittingAlgorithm();
//Report relevant quantities of this and it's splitting planes
	void ReportSplittingAlgorithm();

//Report which splitting object splitted how many particles
	void ReportSplittingSummary(int ncase);

//Do a particle splitting if necessary, very often it will not be needed (when the particle already passed all spliting objects)
	void PotentialParticleSplitting();

//Reset Parameters if necessary
	void PotentialParameterReset();

	virtual void endHistory();

	virtual void EndSimulation(int i_nNumberOfHistories);

	virtual int InitDerivedSplittingAlgorithm(EGS_Input* i_Input);

	virtual std::string WhichOpticalElement();

	virtual int IsInInitialCondition();
	
protected:
	//The last OpticalElement used for path splitting:
	vector<OpticalElement*>::iterator m_pLatestActiveSplittingObject;
	//THe last OpticalElement in the list (which is in reverse order i.e. m_SplittingObject.end() == m_pLatestActiveSplittingObject means that there are no clones around)
	//vector<OpticalElement*>::iterator m_pLastSplittingObject;


	bool m_bCoherentSource;

private:

};

#endif
