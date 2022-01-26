#ifndef _BINARYSOURCEGRATINGFOURIERZPLANE_H_
#define _BINARYSOURCEGRATINGFOURIERZPLANE_H_

#include <fstream>
#include <vector>

#include "egs_rndm.h"

#include "SplittingObject.h"
//#include <vector>

class BinarySourceGratingFourierZPlane : public SplittingObject
{
public:

	BinarySourceGratingFourierZPlane(std::string i_sName, EGS_RandomGenerator* i_pRandomNumberGenerator);
	BinarySourceGratingFourierZPlane(std::string i_sName);
	~BinarySourceGratingFourierZPlane();

	virtual bool IsOnOpticalElement (EGS_Vector* i_pPosition);
	virtual bool ApplyOpticalElementRule();
	virtual void ReportOpticalElement();
	virtual int InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry);
	virtual void PrintSummary();


private:
	EGS_Float m_fPosition;

	EGS_Float m_fHoleSizeInCM;
	EGS_Float m_fPeriodicityInCM;

	EGS_Float m_fTransmissionNorm;

	EGS_Float m_fUpperBoundaryForQx;

	EGS_RandomGenerator* m_pRandomNumberGenerator;

	unsigned int m_nNumberOfDeletedParticles;
};

#endif/*



*/
