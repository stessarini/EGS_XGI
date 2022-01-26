#ifndef _HUYGENSZPLANE_H_
#define _HUYGENSZPLANE_H_

#include "SplittingObject.h"
#include "egs_rndm.h"

class HuygensZPlane : public SplittingObject
{
public:
  HuygensZPlane(std::string i_sName, EGS_RandomGenerator* i_pRandomNumberGenerator);
  virtual ~HuygensZPlane();

  virtual bool ApplyOpticalElementRule();
	virtual void PrintSummary();
	virtual void ReportOpticalElement();
  virtual int InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry);

protected:

  EGS_Float m_fZPosition;
  EGS_Float m_fXmin, m_fXmax;
  EGS_Float m_fYmin, m_fYmax;

  //minimum and maximum angle of split paths w.r.t optical (z) axis
  EGS_Float m_fSplittingAngleMin, m_fSplittingAngleMax;
  EGS_Float m_fRangePolarAngle;

  EGS_RandomGenerator* m_pRandomNumberGenerator;


};
#endif
