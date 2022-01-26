#ifndef _DETECTORZPLANEMEM_H_
#define _DETECTORZPLANEMEM_H_

#include "OpticalElement.h"
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <unordered_map>

struct PixelSignal
{
public:
	EGS_Float RealPart;
	EGS_Float ImaginaryPart;
	int m_nNumberOfPaths;

	PixelSignal()
	{
		RealPart = 0.0;
		ImaginaryPart = 0.0;
		m_nNumberOfPaths = 0;
	}

	PixelSignal(EGS_Float i_fRealPart, EGS_Float i_fImaginaryPart)
	{
		RealPart = i_fRealPart;
		ImaginaryPart = i_fImaginaryPart;
		m_nNumberOfPaths = 1;
	}

	PixelSignal(EGS_Float i_fRealPart, EGS_Float i_fImaginaryPart, int i_nNUmberOfPaths)
	{
		RealPart = i_fRealPart;
		ImaginaryPart = i_fImaginaryPart;
		m_nNumberOfPaths = i_nNUmberOfPaths;
	}

	~PixelSignal()
	{

	}

	void AddSignal(PixelSignal& i_OtherPixelSignal)
	{
		RealPart += i_OtherPixelSignal.RealPart;
		ImaginaryPart += i_OtherPixelSignal.ImaginaryPart;
		m_nNumberOfPaths += i_OtherPixelSignal.m_nNumberOfPaths;
	}

	//Add/score a path amplitude (increase number of paths by 1)
	void AddSignal(EGS_Float& i_fRealPart, EGS_Float i_fImaginaryPart)
	{
		RealPart += i_fRealPart;
		ImaginaryPart += i_fImaginaryPart;
		m_nNumberOfPaths++;
	}

	EGS_Float ReadOut()
	{
		return RealPart * RealPart + ImaginaryPart * ImaginaryPart;
	}

};


class DetectorZPlaneMem : public OpticalElement
{
public:
	DetectorZPlaneMem(std::string i_sName, std::string i_sInputFileName);
	~DetectorZPlaneMem();

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

	std::vector<double>* m_pTotalSignal;
	unsigned int m_nNxTotalSignal;
	unsigned int m_nNyTotalSignal;

	std::unordered_map<unsigned int, PixelSignal> m_mHistorySignal;

	EGS_Vector m_fPosition;
	EGS_Float m_fWidth;
	EGS_Float m_fHeight;

	EGS_Float m_f_TotalSignal_DeltaX; //Bin width
	EGS_Float m_f_Totalsignal_DeltaY;

	bool m_bWirteToBinary;
	bool m_bWriteToPGM;

	//Intput file name
	std::string m_sFileName;
	//path to existing folder. Where to save the results/outputs of the detector
	std::string m_sPath;
	bool m_bAllowedToOverWrite;

	unsigned int m_nNumberOfSecondaryParticles;

	//EGS_Phase_lib_v5_1_1
	bool m_bDoHistoryWiseNormalization;
	EGS_Float m_fTotalIncommingStatWeightOfPrimaryPaths;
	bool m_bDoNumberOfPathsNormalization;

	bool m_bCountNumberOfPaths;
	std::vector<unsigned int>* m_pNumberOfPaths;
};



#endif
