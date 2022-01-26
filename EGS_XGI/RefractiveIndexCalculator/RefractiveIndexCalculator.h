#ifndef _REFRACTIVEINDEXCLCULATOR_H_
#define _REFRACTIVEINDEXCLCULATOR_H_
#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>

#include "egs_base_geometry.h"

class RefractiveIndexCalculator
{
public:
	RefractiveIndexCalculator();
	~RefractiveIndexCalculator();
	void ImportRefractiveIndex(string i_sRefFile, EGS_BaseGeometry* i_pSimulationGeometry);
	EGS_Float GetRefractiveIndexDelta(int i_nMediumIndex, EGS_Float i_fEnergy);
	EGS_Float GetRefractiveIndexDeltaFromWavelength(int i_nMediumIndex, EGS_Float i_fWavelength);


private:
	//The refractive index data that is used to estimate the refractive index with assumed $\lambda^2$-dependence
	std::vector<EGS_Float> m_fRefractiveIndexData;
};


#endif
