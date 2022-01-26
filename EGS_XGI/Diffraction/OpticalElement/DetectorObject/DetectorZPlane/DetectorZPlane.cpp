#include "DetectorZPlane.h"
#include "xgi_global_variables.h"

#include "egs_interface2.h"



DetectorZPlane::DetectorZPlane()
	:OpticalElement(),
	m_pTotalSignal(0),
	m_pHistorySignalRealPart(0),
	m_pHistorySignalImaginaryPart(0),
	m_sPath(""),
	m_nNumberOfSecondaryPhotons(0),
	m_bGenerateSeparateBinaryForSecondarySignal(true),
	m_bDoHistoryWiseNormalization(false),
	m_fTotalIncommingStatWeightOfPrimaryPaths(0.0)
{

}

DetectorZPlane::DetectorZPlane(unsigned int i_nNxTotalSignal, unsigned int i_nNyTotalSignal, EGS_Vector& i_fPosition, EGS_Float i_fWidth, EGS_Float i_fHeight, bool i_bWirteToBinary, bool i_bWriteToPGM, std::string i_sFileName, std::string i_sDetectorName)
	:OpticalElement(),
	m_nNxTotalSignal(i_nNxTotalSignal),
	m_nNyTotalSignal(i_nNyTotalSignal),
	m_fPosition(i_fPosition),
	m_fWidth(i_fWidth),
	m_fHeight(i_fHeight),
	m_bWirteToBinary(i_bWirteToBinary),
	m_bWriteToPGM(i_bWriteToPGM),
	m_sFileName(i_sFileName),
	m_pTotalSignal(0),
	m_pHistorySignalRealPart(0),
	m_pHistorySignalImaginaryPart(0),
	m_pNumberOfPaths(0),
	m_sPath(""),
	m_nNumberOfSecondaryPhotons(0),
	m_bGenerateSeparateBinaryForSecondarySignal(true),
	m_bDoHistoryWiseNormalization(false),
	m_fTotalIncommingStatWeightOfPrimaryPaths(0.0)
{
	m_sName = i_sDetectorName;
	m_f_TotalSignal_DeltaX = m_fWidth/m_nNxTotalSignal;
	m_f_Totalsignal_DeltaY = m_fHeight/m_nNyTotalSignal;
}

DetectorZPlane::DetectorZPlane(std::string i_sName, std::string i_sInputFileName)
	:OpticalElement(),
	m_nNxTotalSignal(1),
	m_nNyTotalSignal(1),
	m_fPosition(),
	m_fWidth(1.0),
	m_fHeight(1.0),
	m_bWirteToBinary(false),
	m_bWriteToPGM(false),
	m_sFileName(i_sInputFileName),
	m_pTotalSignal(0),
	m_pHistorySignalRealPart(0),
	m_pHistorySignalImaginaryPart(0),
	m_pNumberOfPaths(0),
	m_sPath(""),
	m_nNumberOfSecondaryPhotons(0),
	m_bGenerateSeparateBinaryForSecondarySignal(true),
	m_bDoHistoryWiseNormalization(false),
	m_fTotalIncommingStatWeightOfPrimaryPaths(0.0)
{
	m_sName = i_sName;
}

DetectorZPlane::~DetectorZPlane()
{
	if(m_pTotalSignal != 0)
	{
		delete m_pTotalSignal;
	}
	if(m_pHistorySignalRealPart != 0)
	{
		delete m_pHistorySignalRealPart;
	}
	if(m_pHistorySignalImaginaryPart != 0)
	{
		delete m_pHistorySignalImaginaryPart;
	}
	if(m_pNumberOfPaths != 0)
	{
		delete m_pNumberOfPaths;
	}
}

bool DetectorZPlane::IsOnOpticalElement (EGS_Vector* i_pPosition)
{
	if( (abs(i_pPosition->z - m_fPosition.z) < 1e-8) && (m_fPosition.x < i_pPosition->x) && ((m_fPosition.x + m_fWidth) > i_pPosition->x) && (m_fPosition.y < i_pPosition->y) && ((m_fPosition.y + m_fHeight) > i_pPosition->y))
	{
		return true;
	}
	else
	{
		return false;
	}
}



bool DetectorZPlane::ApplyOpticalElementRule()
{
	m_nStackSizeBefore = the_stack->np -1;
	if(abs(the_stack->z[m_nStackSizeBefore] - m_fPosition.z) < 1e-10)
	{
		EGS_Float fXPositionOnDetector = (the_stack->x[m_nStackSizeBefore] - m_fPosition.x);
		EGS_Float fYPositionOnDetector = (the_stack->y[m_nStackSizeBefore] - m_fPosition.y);
		if( (fXPositionOnDetector >= 0.0) && (fXPositionOnDetector < m_fWidth) && (fYPositionOnDetector >= 0.0) && (fYPositionOnDetector < m_fHeight) )
		{
			if(the_stack->iq[m_nStackSizeBefore] == 0 && e_bPrimary[m_nStackSizeBefore] == true)
			{
				m_nNumberOfTimesApplied++;
				EGS_Float fPhase = e_fPhase[m_nStackSizeBefore] * the_stack->E[m_nStackSizeBefore] / ec_fEnergyToWaveLength;
				fPhase = (fPhase - floor(fPhase)) * 2.0 * ec_fPi;//modulo 2 pi
				EGS_Float fnorm = exp(e_fLogNorm[m_nStackSizeBefore]);
				int nDataPointAddres = floor(fXPositionOnDetector/m_f_TotalSignal_DeltaX) + floor(fYPositionOnDetector/m_f_Totalsignal_DeltaY) * m_nNxTotalSignal;
				EGS_Float fSQRT_wt = sqrt(the_stack->wt[m_nStackSizeBefore]);
				(*m_pHistorySignalRealPart)[nDataPointAddres] += fSQRT_wt * fnorm * cos(fPhase);
				(*m_pHistorySignalImaginaryPart)[nDataPointAddres] += fSQRT_wt * fnorm * sin(fPhase);
				(*m_pNumberOfPaths)[nDataPointAddres] += 1;
				if(m_bDoHistoryWiseNormalization == true)
				{
					m_fTotalIncommingStatWeightOfPrimaryPaths += fnorm * fnorm * the_stack->wt[m_nStackSizeBefore];
				}
				return true;
			}
			else if(the_stack->iq[m_nStackSizeBefore] == 0)
			{
				(*m_pTotalSignal)[floor(fXPositionOnDetector/m_f_TotalSignal_DeltaX) + floor(fYPositionOnDetector/m_f_Totalsignal_DeltaY) * m_nNxTotalSignal] += the_stack->wt[m_nStackSizeBefore];
				m_nNumberOfSecondaryPhotons++;
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}

	}
	else
	{
		return false;
	}

}
void DetectorZPlane::PrintSummary()
{
	std::cout << "--------------------------------" << std::endl;
	std::cout << "A detector of type: " << "DetectorZPlane" << std::endl;
	std::cout << "Number scored paths: " << m_nNumberOfTimesApplied << std::endl;
	std::cout << "Number of secondary particles scored: " << m_nNumberOfSecondaryPhotons << std::endl;
	std::cout << "--------------------------------" << std::endl;
}

void DetectorZPlane::ReportOpticalElement()
{
	std::cout << "--------------------------------" << std::endl;
	std::cout << "A detector of type: " << "DetectorZPlane" << std::endl;
	std::cout << "Detector position: (" << m_fPosition.x << ", " << m_fPosition.y << ", " << m_fPosition.z << ")" << std::endl;
	std::cout << "Pixels total signal [Nx, Ny] = [" << m_nNxTotalSignal << ", " << m_nNyTotalSignal << "]" << std::endl;
	std::cout << "Dimensions [x,y] = [" << m_fWidth << ", " << m_fHeight << "]" << std::endl;
	std::cout << "--------------------------------" << std::endl;
}


void DetectorZPlane::WriteOnFile(long int i_nLastCase)
{
	//Setup the binary file

	//try saving it in m_sPath if m_sPath != 0
	bool bBinaryFileIsOpen = false;
	std::string BinaryFileName;
	std::fstream BinarySignalFile;
	std::string AttributeFileName;
	std::fstream AttributeFile;
	if(m_sPath.compare("") != 0)
	{
		BinaryFileName = m_sPath;
		BinaryFileName.append("/");
		BinaryFileName.append(m_sFileName);
		BinaryFileName.append("_");
		BinaryFileName.append(m_sName);
		AttributeFileName = BinaryFileName;
		BinaryFileName.append(".bin");
		BinarySignalFile.open(BinaryFileName, std::ios::binary | std::ios::out | std::ios::trunc);
		bBinaryFileIsOpen = BinarySignalFile.is_open();
		if(!bBinaryFileIsOpen)
		{
			//Failed to open the file.
			//maybe the path doesn't exist or don't have permission
			//try other location
			std::cout << "DetectorZPlane::WriteOnFile: Warning unable to open file: " << BinaryFileName << std::endl;
		}
	}
	if(!bBinaryFileIsOpen)
	{
		BinaryFileName = m_sFileName;
		BinaryFileName.append("/");
		BinaryFileName.append(m_sFileName);
		BinaryFileName.append("_");
		BinaryFileName.append(m_sName);
		AttributeFileName = BinaryFileName;
		BinaryFileName.append(".bin");
		BinarySignalFile.open(BinaryFileName, std::ios::binary | std::ios::out | std::ios::trunc);
		bBinaryFileIsOpen = BinarySignalFile.is_open();
		if(!bBinaryFileIsOpen)
		{
			//Failed to open the file.
			//maybe the path doesn't exist or don't have permission
			//try other location
			std::cout << "DetectorZPlane::WriteOnFile: Warning unable to open file: " << BinaryFileName << std::endl;
		}
	}
	if(!bBinaryFileIsOpen && m_bAllowedToOverWrite)
	{
		BinaryFileName = m_sFileName;
		BinaryFileName.append("_");
		BinaryFileName.append(m_sName);
		AttributeFileName = BinaryFileName;
		BinaryFileName.append(".bin");
		BinarySignalFile.open(BinaryFileName, std::ios::binary | std::ios::out | std::ios::trunc);
		bBinaryFileIsOpen = BinarySignalFile.is_open();
	}
	if(!bBinaryFileIsOpen)
	{
		std::cout << "DetectorZPlane::WriteOnFile: Error unable to open binary file." << std::endl;
	}
	else
	{
		std::cout << "DetectorZPlane::WriteOnFile: " << BinaryFileName << std::endl;
		size_t nTypeSize = sizeof((*m_pTotalSignal)[0]);
		for(vector<double>::iterator PixelIterator = m_pTotalSignal->begin(); PixelIterator != m_pTotalSignal->end(); PixelIterator++)
		{
			BinarySignalFile.write((char*)&(*PixelIterator),nTypeSize);
		}
		BinarySignalFile.close();
	}

	if(bBinaryFileIsOpen)//if opened the binary successfully
	{
		//assuming noone interfered with the file system during the writing process...
		//always want to put the Attribute file at the same location as the binary
		AttributeFileName.append(".att");
		AttributeFile.open(AttributeFileName, std::ios::out | std::ios::trunc);
		if(AttributeFile.is_open())
		{
			AttributeFile << "Corresponding Binary : " << BinaryFileName << std::endl;
			AttributeFile << "Number of histories = " << i_nLastCase << std::endl;
			AttributeFile << "Detector class: " << "DetectorZPlane" << std::endl;
			AttributeFile << "Detector position: (" << m_fPosition.x << ", " << m_fPosition.y << ", " << m_fPosition.z << ")" << std::endl;
			AttributeFile << "Pixels total signal [Nx, Ny] = [" << m_nNxTotalSignal << ", " << m_nNyTotalSignal << "]" << std::endl;
			AttributeFile << "Dimensions [x,y] = [" << m_fWidth << ", " << m_fHeight << "]" << std::endl;
			AttributeFile.close();
		}
	}
}

void DetectorZPlane::WriteToFileNumberOfPaths()
{
	bool bBinaryFileIsOpen = false;
	std::string BinaryFileName;
	std::fstream BinarySignalFile;
	std::string AttributeFileName;
	std::fstream AttributeFile;
	if(m_sPath.compare("") != 0)
	{
		BinaryFileName = m_sPath;
		BinaryFileName.append("/");
		BinaryFileName.append(m_sFileName);
		BinaryFileName.append("_NOP_");
		BinaryFileName.append(m_sName);
		AttributeFileName = BinaryFileName;
		BinaryFileName.append(".bin");
		BinarySignalFile.open(BinaryFileName, std::ios::binary | std::ios::out | std::ios::trunc);
		bBinaryFileIsOpen = BinarySignalFile.is_open();
		if(!bBinaryFileIsOpen)
		{
			//Failed to open the file.
			//maybe the path doesn't exist or don't have permission
			//try other location
			std::cout << "DetectorZPlane::WriteToFileNumberOfPaths: Warning unable to open file: " << BinaryFileName << std::endl;
		}
	}
	if(!bBinaryFileIsOpen)
	{
		BinaryFileName = m_sFileName;
		BinaryFileName.append("/");
		BinaryFileName.append(m_sFileName);
		BinaryFileName.append("_NOP_");
		BinaryFileName.append(m_sName);
		AttributeFileName = BinaryFileName;
		BinaryFileName.append(".bin");
		BinarySignalFile.open(BinaryFileName, std::ios::binary | std::ios::out | std::ios::trunc);
		bBinaryFileIsOpen = BinarySignalFile.is_open();
		if(!bBinaryFileIsOpen)
		{
			//Failed to open the file.
			//maybe the path doesn't exist or don't have permission
			//try other location
			std::cout << "DetectorZPlane::WriteToFileNumberOfPaths: Warning unable to open file: " << BinaryFileName << std::endl;
		}
	}
	if(!bBinaryFileIsOpen && m_bAllowedToOverWrite)
	{
		BinaryFileName = m_sFileName;
		BinaryFileName.append("_NOP_");
		BinaryFileName.append(m_sName);
		AttributeFileName = BinaryFileName;
		BinaryFileName.append(".bin");
		BinarySignalFile.open(BinaryFileName, std::ios::binary | std::ios::out | std::ios::trunc);
		bBinaryFileIsOpen = BinarySignalFile.is_open();
	}
	if(!bBinaryFileIsOpen)
	{
		std::cout << "DetectorZPlane::WriteToFileNumberOfPaths: Error failed to open binary file:" << BinaryFileName << std::endl;
	}
	else
	{
		std::cout << "DetectorZPlane::WriteToFileNumberOfPaths: " << BinaryFileName << std::endl;
		size_t nTypeSize = sizeof((*m_pNumberOfPaths)[0]);
		for(vector<unsigned int>::iterator PixelIterator = m_pNumberOfPaths->begin(); PixelIterator != m_pNumberOfPaths->end(); PixelIterator++)
		{
			BinarySignalFile.write((char*)&(*PixelIterator),nTypeSize);
		}
		if(m_bWirteToBinary == false)
		{
			AttributeFileName.append(".att");
			AttributeFile.open(AttributeFileName, std::ios::out | std::ios::trunc);

			if(AttributeFile.is_open())
			{
				AttributeFile << "Corresponding binary : " << BinaryFileName << std::endl;
				AttributeFile << "Number of scored primary paths = " << m_nNumberOfTimesApplied << std::endl;
				AttributeFile << "Detector class: " << "DetectorZPlaneMem" << std::endl;
				AttributeFile << "Detector position: (" << m_fPosition.x << ", " << m_fPosition.y << ", " << m_fPosition.z << ")" << std::endl;
				AttributeFile << "Pixels total signal [Nx, Ny] = [" << m_nNxTotalSignal << ", " << m_nNyTotalSignal << "]" << std::endl;
				AttributeFile << "Side Lengths [Width, Height] = [" << m_fWidth << ", " << m_fHeight << "]" << std::endl;
				AttributeFile.close();
			}
		}
	}
}

void DetectorZPlane::CreateImage()
{
	EGS_Float fTheMaxValue = 0.0;

	for(std::vector<double>::iterator SignalIterator = m_pTotalSignal->begin(); SignalIterator != m_pTotalSignal->end(); SignalIterator++)
	{
		if(fTheMaxValue < *SignalIterator)
		{
			fTheMaxValue = *SignalIterator;
		}
	}
	if(fTheMaxValue > 0.0)
	{
		//Image header
		ofstream Image;
		std::string name = m_sFileName;
		name.append("/");
		name.append(m_sFileName);
		name.append("_");
		name.append(m_sName);
		name.append(".pgm");
		Image.open(name);

		std::cout << "DetectorZPlane::CreateImage: Write on: " << name << std::endl;

		Image << "P2" << endl << "#comment" << endl << m_nNxTotalSignal << " " << m_nNyTotalSignal << endl << 255 << endl;

		for(unsigned int nYcounter = 0; nYcounter < m_nNyTotalSignal; nYcounter++)
		{
			for (unsigned int nXcounter = 0; nXcounter < m_nNxTotalSignal; nXcounter++)
			{
				Image <<(int)(round(255.0 * (*m_pTotalSignal)[nYcounter * m_nNxTotalSignal + nXcounter] / fTheMaxValue)) << " ";
				if((int)(round(255.0 * (*m_pTotalSignal)[nYcounter * m_nNxTotalSignal + nXcounter] / fTheMaxValue)) == -2147483648)
				{
					std::cout << "pixel: [" << nXcounter << ", " << nYcounter << "] , " << (int)(round(255.0 * (*m_pTotalSignal)[nYcounter * m_nNxTotalSignal + nXcounter] / fTheMaxValue));
				std::cout << "<=> " << (*m_pTotalSignal)[nYcounter * m_nNxTotalSignal + nXcounter] / fTheMaxValue << std::endl;
				}
			}
			Image << endl;
		}
		Image.close();
	}
	else
	{
		std::cout << "WARNING: DetectorZPlane::CreateImage() couldn't create PGM, fTheMaxValue == 0.0" << std::endl;
	}
}


void DetectorZPlane::EndRayTracing()
{
	unsigned nNTotalNumberOfPixels = m_nNxTotalSignal*m_nNyTotalSignal;
	if(m_bGenerateSeparateBinaryForSecondarySignal)
	{
		std::cout << "DetectorZPlane::EndRayTracing print secondary signal" << std::endl;
		std::string copy = m_sFileName;
		m_sFileName.append("_SecondarySignal");
		WriteOnFile(-1);
		m_sFileName = copy;
		//Reset the total signal
		ResetTotalSignal();
	}
	EGS_Float fNormalizationConstant = 1.0;
	if(m_bDoHistoryWiseNormalization == true)
	{
		EGS_Float fTotalSignalGeneratedInThisHistory = 0.0;
		for(unsigned int nPixelCounter = 0; nPixelCounter < nNTotalNumberOfPixels; nPixelCounter++)
		{
			fTotalSignalGeneratedInThisHistory += pow((*m_pHistorySignalRealPart)[nPixelCounter],2) + pow((*m_pHistorySignalImaginaryPart)[nPixelCounter],2);
		}
		fNormalizationConstant = m_fTotalIncommingStatWeightOfPrimaryPaths / fTotalSignalGeneratedInThisHistory;
	}

	for(unsigned int nPixelCounter = 0; nPixelCounter < nNTotalNumberOfPixels; nPixelCounter++)
	{
		(*m_pTotalSignal)[nPixelCounter] += (pow((*m_pHistorySignalRealPart)[nPixelCounter],2) + pow((*m_pHistorySignalImaginaryPart)[nPixelCounter],2)) * fNormalizationConstant;
		//Ressetting the history signal
		(*m_pHistorySignalRealPart)[nPixelCounter] = 0.0;
		(*m_pHistorySignalImaginaryPart)[nPixelCounter] = 0.0;
	}
}


void DetectorZPlane::ScoreSecondaryParticle()
{
	m_nStackSizeBefore = the_stack->np -1;
	if(abs(the_stack->z[m_nStackSizeBefore] - m_fPosition.z < 1e-20))
	{
		EGS_Float fXPositionOnDetector = (the_stack->x[m_nStackSizeBefore] - m_fPosition.x);
		EGS_Float fYPositionOnDetector = (the_stack->y[m_nStackSizeBefore] - m_fPosition.y);
		if( (fXPositionOnDetector >= 0.0) && (fXPositionOnDetector < m_fWidth) && (fYPositionOnDetector >= 0.0) && (fYPositionOnDetector < m_fHeight) )
		{
			m_nNumberOfSecondaryPhotons++;
			(*m_pTotalSignal)[floor(fXPositionOnDetector/m_f_TotalSignal_DeltaX) + floor(fYPositionOnDetector/m_f_Totalsignal_DeltaY) * m_nNxTotalSignal] += the_stack->wt[m_nStackSizeBefore];
		}
	}
}

int DetectorZPlane::InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry)
{
	std::vector<int> nTotalSignalPixel;
	std::vector<EGS_Float> fPosition;
	std::vector<EGS_Float> fHeight;
	std::vector<EGS_Float> fWidth;
	std::vector<int> nForBoolWirteToBinary;
	std::vector<int> nForBoolWirteToPGM;
	std::string sCountNumberOfPaths;
	std::string sPath;
	std::string sAllowedToOverWrite;
	std::string sGenerateSecondarySignalFile;
	std::string sDoHistoryWiseNormalization;

	int errTotSig = i_sInput->getInput("total signal pixels",nTotalSignalPixel);
	int errPos = i_sInput->getInput("position",fPosition);
	int errHeight = i_sInput->getInput("height",fHeight);
	int errWidth = i_sInput->getInput("width",fWidth);
	int errBinary = i_sInput->getInput("write to binary",nForBoolWirteToBinary);
	int errPGM = i_sInput->getInput("write to pgm",nForBoolWirteToPGM);
	int errNPath = i_sInput->getInput("count number of paths", sCountNumberOfPaths);
	int errsPath = i_sInput->getInput("save at", sPath);
	int errsOW = i_sInput->getInput("allowed to overwrite", sAllowedToOverWrite);
	int errSec = i_sInput->getInput("generate secondary signal file", sGenerateSecondarySignalFile);
	int errHistNorm = i_sInput->getInput("history wise normalization", sDoHistoryWiseNormalization);

	if(errNPath)
	{
		std::cout << "DetectorZPlane::InitOpticalElement: Error reading 'count number of paths'-key. Use default - turn on." << std::endl;
		m_bCountNumberOfPaths = true;
	}
	else
	{
		if(sCountNumberOfPaths.compare("yes") == 0)
		{
			m_bCountNumberOfPaths = true;
		}
		else
		{
			if(sCountNumberOfPaths.compare("no") != 0)
			{
				std::cout << "DetectorZPlane::InitOpticalElement: WARNING key for flag 'count number of paths' not understood. Use: 'yes' or 'no'. Use default: output number of paths." << std::endl;
				m_bCountNumberOfPaths = true;
			}
			else
			{
				m_bCountNumberOfPaths = false;
			}
		}
	}
	if(errsPath)
	{
		m_sPath = "";
		std::cout << "DetectorZPlane::InitOpticalElement: Error reading 'save at'. Save results in './m_sFileName/' or '.'." << std::endl;
	}
	else
	{
		m_sPath = sPath;
	}
	if(errsOW)
	{
		std::cout << "******************************************************************************************" << std::endl;
		std::cout << "DetectorZPlane::InitOpticalElement: Error reading 'allowed to overwrite'. DON'T OVERWRITE." << std::endl;
		std::cout << "******************************************************************************************" << std::endl;
		m_bAllowedToOverWrite = false;
	}
	else
	{
		if(sAllowedToOverWrite.compare("yes") == 0)
		{
			m_bAllowedToOverWrite = true;
		}
		else if(sAllowedToOverWrite.compare("no") == 0)
		{
			m_bAllowedToOverWrite = false;
		}
		else
		{
			std::cout << "******************************************************************************************" << std::endl;
			std::cout << "DetectorZPlane::InitOpticalElement: 'allowed to overwrite'-key unknown. DON'T OVERWRITE." << std::endl;
			std::cout << "******************************************************************************************" << std::endl;
			m_bAllowedToOverWrite = false;
		}
	}
	if(errHistNorm)
	{
		std::cout << "DetectorZPlaneMem::InitOpticalElement: Warning error reading flag 'history wise normalization'. Use default: Disable history wise normalization." << std::endl;
		m_bDoHistoryWiseNormalization = false;
	}
	else
	{
		if(sDoHistoryWiseNormalization.compare("yes") == 0)
		{
			m_bDoHistoryWiseNormalization = true;
		}
		else if(sDoHistoryWiseNormalization.compare("no") == 0)
		{
			m_bDoHistoryWiseNormalization = false;
		}
		else
		{
			std::cout << "DetectorZPlaneMem::InitOpticalElement: Warning unknown value for 'history wise normalization'. Use default: Disable history wise normalization." << std::endl;
			m_bDoHistoryWiseNormalization = false;
		}
	}
	if(errTotSig || nTotalSignalPixel.size() != 2)
	{
		std::cout << "DetectorZPlane::InitOpticalElement: Error reading number of pixels of total signal." << std::endl;
		return 0;
	}
	else if(errPos || fPosition.size() != 3)
	{
		std::cout << "DetectorZPlane::InitOpticalElement: Error reading position." << std::endl;
		return 0;
	}
	else if(errHeight || fHeight.size() != 1)
	{
		std::cout << "DetectorZPlane::InitOpticalElement: Error reading height." << std::endl;
		return 0;
	}
	else if(errWidth || fWidth.size() != 1)
	{
		std::cout << "DetectorZPlane::InitOpticalElement: Error reading height." << std::endl;
		return 0;
	}
	else if(errBinary || nForBoolWirteToBinary.size() != 1)
	{
		std::cout << "DetectorZPlane::InitOpticalElement: Error reading flag for binary output." << std::endl;
		return 0;
	}
	else if(errPGM || nForBoolWirteToPGM.size() != 1)
	{
		std::cout << "DetectorZPlane::InitOpticalElement: Error reading flag for PGM output." << std::endl;
		return 0;
	}
	else
	{
		if(errSec)
		{
			std::cout << "DetectorZPlane::InitOpticalElement: Unable to read 'generate secondary signal file'. Disable generation of secondary signal file." << std::endl;
			m_bGenerateSeparateBinaryForSecondarySignal = false;
		}
		else if(sGenerateSecondarySignalFile == "yes")
		{
			m_bGenerateSeparateBinaryForSecondarySignal = true;
		}
		else if(sGenerateSecondarySignalFile == "no")
		{
			m_bGenerateSeparateBinaryForSecondarySignal = false;
		}
		else
		{
			std::cout << "DetectorZPlane::InitOpticalElement: Unknown option for 'generate secondary signal file'. Disable generation of secondary signal file." << std::endl;
			m_bGenerateSeparateBinaryForSecondarySignal = false;
		}
		m_nNxTotalSignal = nTotalSignalPixel[0];
		m_nNyTotalSignal = nTotalSignalPixel[1];
		m_fPosition.x = fPosition[0];
		m_fPosition.y = fPosition[1];
		m_fPosition.z = fPosition[2];
		m_fWidth = fWidth[0];
		m_fHeight = fHeight[0];

		m_f_TotalSignal_DeltaX = m_fWidth/m_nNxTotalSignal;
		m_f_Totalsignal_DeltaY = m_fHeight/m_nNyTotalSignal;

		m_pTotalSignal = new std::vector<double> (m_nNxTotalSignal * m_nNyTotalSignal,0.0);
		m_pHistorySignalRealPart = new std::vector<double> (m_nNxTotalSignal * m_nNyTotalSignal,0.0);
		m_pHistorySignalImaginaryPart = new std::vector<double> (m_nNxTotalSignal * m_nNyTotalSignal,0.0);
		m_pNumberOfPaths = new std::vector<unsigned int> (m_nNxTotalSignal * m_nNyTotalSignal, 0);

		if(nForBoolWirteToBinary[0] == 0)
		{
			m_bWirteToBinary = false;
		}
		else
		{
			m_bWirteToBinary = true;
		}
		if(nForBoolWirteToPGM[0] == 0)
		{
			m_bWriteToPGM = false;
		}
		else
		{
			m_bWriteToPGM = true;
		}
	}
	return 1;
}

void DetectorZPlane::EndSimulation(int i_nNumberOfHistories)
{
	if(m_bWirteToBinary == true)
	{
		WriteOnFile(i_nNumberOfHistories);
	}
	if(m_bWriteToPGM == true)
	{
		CreateImage();
	}
	if(m_bCountNumberOfPaths == true)
	{
		WriteToFileNumberOfPaths();
	}
	std::cout << "DetectorZPlane::EndSimulation(): scored primary paths: " << m_nNumberOfTimesApplied << std::endl;
	std::cout << "DetectorZPlane::EndSimulation(): scored secondary photons " << m_nNumberOfSecondaryPhotons << std::endl;
}


bool DetectorZPlane::PermissionToResetPointer(int i_nStackSize)
{
	return true;
}


void DetectorZPlane::ResetTotalSignal()
{
	for(std::vector<double>::iterator TotalSignalIterator = m_pTotalSignal->begin(); TotalSignalIterator != m_pTotalSignal->end(); TotalSignalIterator++)
	{
		*TotalSignalIterator = 0.0;
	}
}
