#ifndef _DETECTORZPLANE_H_
#define _DETECTORZPLANE_H_

#include "OpticalElement.h"
#include <vector>
#include <string>
#include <fstream>


class DetectorZPlane : public OpticalElement
{
public:
	DetectorZPlane();

	DetectorZPlane(unsigned int i_nNxTotalSignal, unsigned int i_nNyTotalSignal, EGS_Vector& i_fPosition, EGS_Float i_fWidth, EGS_Float i_fHeight, bool i_bWirteToBinary, bool i_bWriteToPGM, std::string i_sFileName, std::string i_sDetectorName);

	DetectorZPlane(std::string i_sName, std::string i_sInputFileName);

	~DetectorZPlane();

	virtual bool IsOnOpticalElement (EGS_Vector* i_pPosition);
	virtual bool ApplyOpticalElementRule();
	virtual void PrintSummary();
	virtual void ReportOpticalElement();

	virtual void EndRayTracing();

	void WriteOnFile(long int i_nLastCase);
	void CreateImage();
	void WriteToFileNumberOfPaths();
	virtual void ScoreSecondaryParticle();
	virtual bool PermissionToResetPointer(int i_nStackSize);
	virtual int InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry);

	virtual void EndSimulation(int i_nNumberOfHistories);
protected:

	void ResetTotalSignal();

	std::vector<double>* m_pTotalSignal;
	std::vector<double>* m_pHistorySignalRealPart;
	std::vector<double>* m_pHistorySignalImaginaryPart;

	unsigned int m_nNxTotalSignal;
	unsigned int m_nNyTotalSignal;

	//unsigned int m_nNxHistorySignal;
	//unsigned int m_nNyHistorySignal;

	EGS_Vector m_fPosition;
	EGS_Float m_fWidth;
	EGS_Float m_fHeight;

	EGS_Float m_f_TotalSignal_DeltaX; //Bin width
	EGS_Float m_f_Totalsignal_DeltaY;

	//EGS_Float m_f_HistorySignal_DeltaX;
	//EGS_Float m_f_HistorySignal_DeltaY;

	bool m_bWirteToBinary;
	bool m_bWriteToPGM;

	std::string m_sFileName;

	bool m_bCountNumberOfPaths;
	std::vector<unsigned int>* m_pNumberOfPaths;

	std::string m_sPath;
	bool m_bAllowedToOverWrite;

	unsigned int m_nNumberOfSecondaryPhotons;

	bool m_bGenerateSeparateBinaryForSecondarySignal;
	bool m_bDoHistoryWiseNormalization;
	EGS_Float m_fTotalIncommingStatWeightOfPrimaryPaths;
};

#endif
