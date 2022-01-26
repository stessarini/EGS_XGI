#ifndef _SPLITTINGOBJECT_H_
#define _SPLITTINGOBJECT_H_


#include <string>

#include "OpticalElement.h"


class SplittingObject : public OpticalElement
{
public:


	SplittingObject(std::string i_sName, std::string i_sType);
	virtual ~SplittingObject();

	virtual bool IsOnOpticalElement (EGS_Vector* i_pPosition);
	virtual bool ApplyOpticalElementRule();
	virtual void PrintSummary();
	virtual void ReportOpticalElement();

	virtual bool PermissionToResetPointer(int i_nStackSize);

	std::string GetName();
	std::string GetType();

	int GetNumberOfSplittingEvents();
	int GetNumberOfGeneratedPaths();


protected:

	int m_nMinimumStackSizeAfter;
  std::string m_sType;
  unsigned int m_nNumberOfGeneratedPaths;
  int m_nSplittingNumber;
};

#endif
