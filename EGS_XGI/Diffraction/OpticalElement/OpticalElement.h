#ifndef _OPTICALELEMENT_H_
#define _OPTICALELEMENT_H_

#include "egs_vector.h"
#include "egs_input.h"

class EGS_Interpolator;
class EGS_BaseGeometry;
class OpticalElement
{
public:
	OpticalElement();
	virtual ~OpticalElement();
	virtual bool IsOnOpticalElement (EGS_Vector* i_pPosition);
	virtual bool ApplyOpticalElementRule();
	virtual void PrintSummary();
	virtual void ReportOpticalElement();
	//the following gives the opportunity to precalculate parameters that are constant untill EndRayTracing is called and the energy of primary photons has changed
	virtual void StartRayTracing(EGS_Float& i_fEnergy);
	//A function to call when all paths of a photon are transported. -> Add detector signals, reset phase gratings
	virtual void EndRayTracing();
	virtual void ScoreSecondaryParticle();
	virtual int InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry);

	//Give the SplittingAlgorithm the permission to reset its pointer to the previous OpticalElement if all paths created by this are transported
	virtual bool PermissionToResetPointer(int i_nStackSize);
	EGS_Float GetNumberOfTimesAppliedOpticalElemetRule();
	virtual void EndSimulation(int i_nNumberOfHistories);

	std::string GetName /*const*/();

	double GetPhaseOfComplexNumber(EGS_Float i_fRealPart, EGS_Float i_fImaginaryPart);

protected:
	unsigned int m_nNumberOfTimesApplied;

	std::string m_sName;

	int m_nStackSizeBefore;
private:

};

#endif
