#include "DetectorZPlaneMem.h"

#include "xgi_global_variables.h"


#include "egs_interface2.h"

DetectorZPlaneMem::DetectorZPlaneMem(std::string i_sName, std::string i_sInputFileName)
	:OpticalElement(),
	m_sFileName(i_sInputFileName),
	m_pTotalSignal(0),
	m_mHistorySignal(),
	m_nNumberOfSecondaryParticles(0),
	m_bDoHistoryWiseNormalization(false),
	m_fTotalIncommingStatWeightOfPrimaryPaths(0.0),
	m_bDoNumberOfPathsNormalization(false),
	m_bCountNumberOfPaths(false),
	m_pNumberOfPaths(0),
	m_sPath(""),
	m_bAllowedToOverWrite(false)
{
	m_sName = i_sName;
}

DetectorZPlaneMem::~DetectorZPlaneMem()
{
	if(m_pTotalSignal != 0)
	{
		delete m_pTotalSignal;
	}
	if(m_pNumberOfPaths != 0)
	{
		delete m_pNumberOfPaths;
	}
}

bool DetectorZPlaneMem::IsOnOpticalElement (EGS_Vector* i_pPosition)
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

bool DetectorZPlaneMem::ApplyOpticalElementRule()
{
	m_nStackSizeBefore = the_stack->np -1;

	if(abs(the_stack->z[m_nStackSizeBefore] - m_fPosition.z) < 1e-20)
	{


		EGS_Float fXPositionOnDetector = (the_stack->x[m_nStackSizeBefore] - m_fPosition.x);
		EGS_Float fYPositionOnDetector = (the_stack->y[m_nStackSizeBefore] - m_fPosition.y);

		if( (fXPositionOnDetector >= 0.0) && (fXPositionOnDetector < m_fWidth) && (fYPositionOnDetector >= 0.0) && (fYPositionOnDetector < m_fHeight) )
		{
			if(the_stack->iq[m_nStackSizeBefore] == 0 && e_bPrimary[m_nStackSizeBefore] == true)
			{
				m_nNumberOfTimesApplied++;

				EGS_Float fPhase = e_fPhase[m_nStackSizeBefore] * the_stack->E[m_nStackSizeBefore] /ec_fEnergyToWaveLength;
				fPhase = (fPhase - floor(fPhase)) * 2.0 * ec_fPi;//modulo 2 pi
				EGS_Float fnorm = exp(e_fLogNorm[m_nStackSizeBefore]);

				//the signal
				EGS_Float fRealPart = sqrt(the_stack->wt[m_nStackSizeBefore]) * fnorm * cos(fPhase);
				EGS_Float fImaginaryPart = sqrt(the_stack->wt[m_nStackSizeBefore]) * fnorm * sin(fPhase);

				unsigned int nPixelAddress = floor(fXPositionOnDetector/m_f_TotalSignal_DeltaX) + floor(fYPositionOnDetector/m_f_Totalsignal_DeltaY) * m_nNxTotalSignal;
				PixelSignal PathSignal(fRealPart, fImaginaryPart);
				//It appears that one can not use emplace() with gcc 447
				//pair<map<unsigned int, PixelSignal>::iterator, bool> pPixelToWriteOn = m_mHistorySignal.emplace(nPixelAddress, a);

				//Add a new pixel to the history signal if that pixel didn't already receive some signal from the current history...
				pair<std::unordered_map<unsigned int, PixelSignal>::iterator, bool> pPixelToWriteOn = m_mHistorySignal.insert(std::pair<unsigned int, PixelSignal>(nPixelAddress, PathSignal));

				if(pPixelToWriteOn.second == false)
				{
					(pPixelToWriteOn.first->second).AddSignal(PathSignal);
				}
				if(m_bDoHistoryWiseNormalization == true)
				{
					m_fTotalIncommingStatWeightOfPrimaryPaths += fnorm * fnorm * the_stack->wt[m_nStackSizeBefore];
				}
				return true;
			}
			else if(the_stack->iq[m_nStackSizeBefore] == 0)
			{
				m_nNumberOfSecondaryParticles++;
				(*m_pTotalSignal)[floor(fXPositionOnDetector/m_f_TotalSignal_DeltaX) + floor(fYPositionOnDetector/m_f_Totalsignal_DeltaY) * m_nNxTotalSignal] += the_stack->wt[m_nStackSizeBefore];
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

void DetectorZPlaneMem::PrintSummary()
{
	std::cout << "--------------------------------" << std::endl;
	std::cout << "A detector of type: " << "DetectorZPlaneMem" << std::endl;
	std::cout << "Number scored paths: " << m_nNumberOfTimesApplied << std::endl;
	std::cout << "Number of secondary particles: " << m_nNumberOfSecondaryParticles << std::endl;
	std::cout << "--------------------------------" << std::endl;
}

void DetectorZPlaneMem::ReportOpticalElement()
{
	std::cout << "--------------------------------" << std::endl;
	std::cout << "A detector of type: " << "DetectorZPlaneMem" << std::endl;
	std::cout << "Detector position: (" << m_fPosition.x << ", " << m_fPosition.y << ", " << m_fPosition.z << ")" << std::endl;
	std::cout << "Pixels total signal [Nx, Ny] = [" << m_nNxTotalSignal << ", " << m_nNyTotalSignal << "]" << std::endl;
	std::cout << "Dimensions [x,y] = [" << m_fWidth << ", " << m_fHeight << "]" << std::endl;
	std::cout << "Has history wise signal normalization: ";
	if(m_bDoHistoryWiseNormalization == true)
	{
		std::cout << "yes" << std::endl;
	}
	else
	{
		std::cout << "no" << std::endl;
	}
	std::cout << "Number of paths normalization: ";
	if(m_bDoNumberOfPathsNormalization == true)
	{
		std::cout << "on" << std::endl;
	}
	else
	{
		std::cout << "off" << std::endl;
	}
	std::cout << "--------------------------------" << std::endl;
}

void DetectorZPlaneMem::WriteOnFile(long int i_nLastCase)
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
			std::cout << "DetectorZPlaneMem::WriteOnFile: Warning unable to open file: " << BinaryFileName << std::endl;
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
			std::cout << "DetectorZPlaneMem::WriteOnFile: Warning unable to open file: " << BinaryFileName << std::endl;
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
		std::cout << "DetectorZPlaneMem::WriteOnFile: Error unable to open binary file." << std::endl;
	}
	else
	{
		std::cout << "DetectorZPlaneMem::WriteOnFile: " << BinaryFileName << std::endl;
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
			AttributeFile << "Number of scored primary paths = " << m_nNumberOfTimesApplied << std::endl;
			AttributeFile << "Detector class: " << "DetectorZPlaneMem" << std::endl;
			AttributeFile << "Detector position: (" << m_fPosition.x << ", " << m_fPosition.y << ", " << m_fPosition.z << ")" << std::endl;
			AttributeFile << "Pixels total signal [Nx, Ny] = [" << m_nNxTotalSignal << ", " << m_nNyTotalSignal << "]" << std::endl;
			AttributeFile << "Dimensions [x,y] = [" << m_fWidth << ", " << m_fHeight << "]" << std::endl;
			AttributeFile << "history wise normalization: ";
			if(m_bDoHistoryWiseNormalization == true)
			{
				AttributeFile << "on" << std::endl;
			}
			else
			{
				AttributeFile << "off" << std::endl;
			}
			AttributeFile << "number of paths normalization: ";
			if(m_bDoNumberOfPathsNormalization == true)
			{
				AttributeFile << "on" << std::endl;
			}
			else
			{
				AttributeFile << "off" << std::endl;
			}
			AttributeFile.close();
		}
	}
}

void DetectorZPlaneMem::WriteToFileNumberOfPaths()
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
			std::cout << "DetectorZPlaneMem::WriteToFileNumberOfPaths: Warning unable to open file: " << BinaryFileName << std::endl;
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
			std::cout << "DetectorZPlaneMem::WriteToFileNumberOfPaths: Warning unable to open file: " << BinaryFileName << std::endl;
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
		std::cout << "DetectorZPlaneMem::WriteToFileNumberOfPaths: Error failed to open binary file:" << BinaryFileName << std::endl;
	}
	else
	{
		std::cout << "DetectorZPlaneMem::WriteToFileNumberOfPaths: " << BinaryFileName << std::endl;
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

void DetectorZPlaneMem::CreateImage()
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
		//std::cout << "DetectorZPlaneMem::CreateImage::fTheMaxValue: " << fTheMaxValue << std::endl;
		ofstream Image;

		bool bImageFileIsOpen = false;
		std::string name;
		if(m_sPath.compare("") != 0)
		{
			name = m_sPath;
			name.append("/");
			name.append(m_sFileName);
			name.append("_");
			name.append(m_sName);
			name.append(".pgm");
			Image.open(name);
			bImageFileIsOpen = Image.is_open();
			if(!bImageFileIsOpen)
			{
				//Failed to open the file.
				//maybe the path doesn't exist or don't have permission
				//try other location
				std::cout << "DetectorZPlaneMem::CreateImage: Warning unable to open file: " << name << std::endl;
			}
		}
		if(!bImageFileIsOpen)
		{
			name = m_sFileName;
			name.append("/");
			name.append(m_sFileName);
			name.append("_");
			name.append(m_sName);
			name.append(".pgm");
			Image.open(name);
			bImageFileIsOpen = Image.is_open();
			if(!bImageFileIsOpen)
			{
				//Failed to open the file.
				//maybe the path doesn't exist or don't have permission
				//try other location
				std::cout << "DetectorZPlaneMem::CreateImage: Warning unable to open file: " << name << std::endl;
			}
		}
		if(!bImageFileIsOpen && m_bAllowedToOverWrite)
		{
			name = m_sFileName;
			name.append("_");
			name.append(m_sName);
			name.append(".pgm");
			Image.open(name);
			bImageFileIsOpen = Image.is_open();
		}
		if(!bImageFileIsOpen)
		{
			std::cout << "DetectorZPlaneMem::CreateImage: Error failed to open file." << std::endl;
		}
		else
		{
			std::cout << "DetectorZPlaneMem::CreateImage: Write on: " << name << std::endl;
			Image << "P2" << endl << "#comment" << endl << m_nNxTotalSignal << " " << m_nNyTotalSignal << endl << 255 << endl;
			for(unsigned int nYcounter = 0; nYcounter < m_nNyTotalSignal; nYcounter++)
			{
				for (unsigned int nXcounter = 0; nXcounter < m_nNxTotalSignal; nXcounter++)
				{
					Image <<(int)(round(255.0 * (*m_pTotalSignal)[nYcounter * m_nNxTotalSignal + nXcounter] / fTheMaxValue)) << " ";
				}
				Image << endl;
			}
			Image.close();
		}
	}
	else
	{
		std::cout << "WARNING: DetectorZPlaneMem::CreateImage() couldn't create PGM, fTheMaxValue == 0.0" << std::endl;
	}
}

void DetectorZPlaneMem::EndRayTracing()
{
	EGS_Float fNormalizationConstant = 1.0;
	//The summing over paths is only allowed to shift intensity to another place. The over all intensity on the detector has to remain the same
	if(m_bDoHistoryWiseNormalization == true && m_bDoNumberOfPathsNormalization == false)
	{
		EGS_Float fTotalSignalGeneratedInThisHistory = 0.0;
		for(std::unordered_map<unsigned int, PixelSignal>::iterator PixelIterator = m_mHistorySignal.begin(); PixelIterator != m_mHistorySignal.end(); PixelIterator++)
		{
			fTotalSignalGeneratedInThisHistory += PixelIterator->second.ReadOut();
		}
		if(fTotalSignalGeneratedInThisHistory > 0.0)
		{
			fNormalizationConstant = m_fTotalIncommingStatWeightOfPrimaryPaths / fTotalSignalGeneratedInThisHistory;
		}
		for(std::unordered_map<unsigned int, PixelSignal>::iterator PixelIterator = m_mHistorySignal.begin(); PixelIterator != m_mHistorySignal.end(); PixelIterator++)
		{
			(*m_pTotalSignal)[PixelIterator->first] += PixelIterator->second.ReadOut() * fNormalizationConstant;
		}
	}
	else if(m_bDoHistoryWiseNormalization == true && m_bDoNumberOfPathsNormalization == true)
	{
		EGS_Float fTotalSignalGeneratedInThisHistory = 0.0;
		for(std::unordered_map<unsigned int, PixelSignal>::iterator PixelIterator = m_mHistorySignal.begin(); PixelIterator != m_mHistorySignal.end(); PixelIterator++)
		{
			fTotalSignalGeneratedInThisHistory += PixelIterator->second.ReadOut() / (PixelIterator->second.m_nNumberOfPaths * PixelIterator->second.m_nNumberOfPaths);
		}
		if(fTotalSignalGeneratedInThisHistory > 0.0)
		{
			fNormalizationConstant = m_fTotalIncommingStatWeightOfPrimaryPaths / fTotalSignalGeneratedInThisHistory;
		}
		for(std::unordered_map<unsigned int, PixelSignal>::iterator PixelIterator = m_mHistorySignal.begin(); PixelIterator != m_mHistorySignal.end(); PixelIterator++)
		{
			(*m_pTotalSignal)[PixelIterator->first] += PixelIterator->second.ReadOut() * fNormalizationConstant / (PixelIterator->second.m_nNumberOfPaths * PixelIterator->second.m_nNumberOfPaths);
		}
	}
	else if(m_bDoNumberOfPathsNormalization == true)
	{
		for(std::unordered_map<unsigned int, PixelSignal>::iterator PixelIterator = m_mHistorySignal.begin(); PixelIterator != m_mHistorySignal.end(); PixelIterator++)
		{
			(*m_pTotalSignal)[PixelIterator->first] += PixelIterator->second.ReadOut() * fNormalizationConstant / (PixelIterator->second.m_nNumberOfPaths * PixelIterator->second.m_nNumberOfPaths);
		}
	}
	else
	{
		for(std::unordered_map<unsigned int, PixelSignal>::iterator PixelIterator = m_mHistorySignal.begin(); PixelIterator != m_mHistorySignal.end(); PixelIterator++)
		{
			(*m_pTotalSignal)[PixelIterator->first] += PixelIterator->second.ReadOut();
		}
	}
	if(m_bCountNumberOfPaths == true)
	{
		for(std::unordered_map<unsigned int, PixelSignal>::iterator PixelIterator = m_mHistorySignal.begin(); PixelIterator != m_mHistorySignal.end(); PixelIterator++)
		{
			(*m_pNumberOfPaths)[PixelIterator->first] += PixelIterator->second.m_nNumberOfPaths;
		}
	}
	m_fTotalIncommingStatWeightOfPrimaryPaths = 0.0;
	m_mHistorySignal.clear();
}

void DetectorZPlaneMem::ScoreSecondaryParticle()
{
	m_nStackSizeBefore = the_stack->np -1;
	if(abs(the_stack->z[m_nStackSizeBefore] - m_fPosition.z < 1e-20))
	{
		EGS_Float fXPositionOnDetector = (the_stack->x[m_nStackSizeBefore] - m_fPosition.x);
		EGS_Float fYPositionOnDetector = (the_stack->y[m_nStackSizeBefore] - m_fPosition.y);
		if( (fXPositionOnDetector >= 0.0) && (fXPositionOnDetector < m_fWidth) && (fYPositionOnDetector >= 0.0) && (fYPositionOnDetector < m_fHeight) )
		{
			m_nNumberOfSecondaryParticles++;
			(*m_pTotalSignal)[floor(fXPositionOnDetector/m_f_TotalSignal_DeltaX) + floor(fYPositionOnDetector/m_f_Totalsignal_DeltaY) * m_nNxTotalSignal] += the_stack->wt[m_nStackSizeBefore];
		}
	}
}

int DetectorZPlaneMem::InitOpticalElement(EGS_Input* i_sInput, EGS_BaseGeometry* i_pGeometry)
{
	std::vector<int> nTotalSignalPixel;
	std::vector<EGS_Float> fPosition;
	std::vector<EGS_Float> fHeight;
	std::vector<EGS_Float> fWidth;
	std::vector<int> nForBoolWirteToBinary;
	std::vector<int> nForBoolWirteToPGM;
	std::string sDoHistoryWiseNormalization;
	std::string sDoNumberOfPathsNormalization;
	std::string sCountNumberOfPaths;
	std::string sPath;
	std::string sAllowedToOverWrite;

	int errTotSig = i_sInput->getInput("total signal pixels",nTotalSignalPixel);
	int errPos = i_sInput->getInput("position",fPosition);
	int errHeight = i_sInput->getInput("height",fHeight);
	int errWidth = i_sInput->getInput("width",fWidth);
	int errBinary = i_sInput->getInput("write to binary",nForBoolWirteToBinary);
	int errPGM = i_sInput->getInput("write to pgm",nForBoolWirteToPGM);
	int errHistNorm = i_sInput->getInput("history wise normalization", sDoHistoryWiseNormalization);
	int errNumNorm = i_sInput->getInput("number of paths normalization", sDoNumberOfPathsNormalization);
	int errNPath = i_sInput->getInput("count paths", sCountNumberOfPaths);
	int errsPath = i_sInput->getInput("save at", sPath);
	int errsOW = i_sInput->getInput("allowed to overwrite", sAllowedToOverWrite);

	if(errTotSig || nTotalSignalPixel.size() != 2)
	{
		std::cout << "DetectorZPlaneMem::InitOpticalElement: Error reading number of pixels of total signal." << std::endl;
		return 0;
	}
	else if(errPos || fPosition.size() != 3)
	{
		std::cout << "DetectorZPlaneMem::InitOpticalElement: Error reading position." << std::endl;
		return 0;
	}
	else if(errHeight || fHeight.size() != 1)
	{
		std::cout << "DetectorZPlaneMem::InitOpticalElement: Error reading height." << std::endl;
		return 0;
	}
	else if(errWidth || fWidth.size() != 1)
	{
		std::cout << "DetectorZPlaneMem::InitOpticalElement: Error reading height." << std::endl;
		return 0;
	}
	else if(errBinary || nForBoolWirteToBinary.size() != 1)
	{
		std::cout << "DetectorZPlaneMem::InitOpticalElement: Error reading flag for binary output." << std::endl;
		return 0;
	}
	else if(errPGM || nForBoolWirteToPGM.size() != 1)
	{
		std::cout << "DetectorZPlaneMem::InitOpticalElement: Error reading flag for PGM output." << std::endl;
		return 0;
	}
	else
	{
		if (errHistNorm)
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
			else
			{
				m_bDoHistoryWiseNormalization = false;
			}
		}
		if(errNumNorm)
		{
			std::cout << "DetectorZPlaneMem::InitOpticalElement: Warning error reading flag 'number of paths normalization'. Use default: Disable number of paths normalization." << std::endl;
			m_bDoNumberOfPathsNormalization = false;
		}
		else
		{
			if(sDoNumberOfPathsNormalization.compare("yes") == 0)
			{
				m_bDoNumberOfPathsNormalization = true;
			}
			else
			{
				if(sDoNumberOfPathsNormalization.compare("no") != 0)
				{
					std::cout << "DetectorZPlaneMem::InitOpticalElement: WARNING key for flag 'number of paths normalization' not understood. Use: 'yes' or 'no'. Use default: Disable number of paths normalization." << std::endl;
				}
				m_bDoNumberOfPathsNormalization = false;
			}
		}
		if(errNPath)
		{
			std::cout << "DetectorZPlaneMem::InitOpticalElement: Warning 1: error reading flag 'count paths'. Use default: Do not count total number of paths." << std::endl;
			m_bCountNumberOfPaths == false;
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
					std::cout << "DetectorZPlaneMem::InitOpticalElement: WARNING 2: key for flag 'count paths' not understood. Use: 'yes' or 'no'. Use default: Do not count total number of paths." << std::endl;
				}
				m_bCountNumberOfPaths = false;
			}
		}
		if(errsPath)
		{
			m_sPath = "";
			std::cout << "DetectorZPlaneMem::InitOpticalElement: Error reading 'save at'. Save results in './m_sFileName/' or '.'." << std::endl;
		}
		else
		{
			m_sPath = sPath;
		}
		if(errsOW)
		{
			std::cout << "******************************************************************************************" << std::endl;
			std::cout << "DetectorZPlaneMem::InitOpticalElement: Error reading 'allowed to overwrite'. DON'T OVERWRITE." << std::endl;
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
				std::cout << "DetectorZPlaneMem::InitOpticalElement: 'allowed to overwrite'-key unknown. DON'T OVERWRITE." << std::endl;
				std::cout << "******************************************************************************************" << std::endl;
				m_bAllowedToOverWrite = false;
			}
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


		//Setup the signal arrays with 0 entries
		m_pTotalSignal = new std::vector<double> (m_nNxTotalSignal*m_nNyTotalSignal, 0.0);
		if(m_bCountNumberOfPaths == true)
		{
			m_pNumberOfPaths = new std::vector<unsigned int> (m_nNxTotalSignal*m_nNyTotalSignal, 0);
		}


		if(nForBoolWirteToBinary[0] == 0)//c++ convention
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


void DetectorZPlaneMem::EndSimulation(int i_nNumberOfHistories)
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
	//std::cout << "DetectorZPlaneMem::EndSimulation(): scored primary paths: " << m_nNumberOfTimesApplied << std::endl;
	//std::cout << "DetectorZPlaneMem::EndSimulation(): scored secondary paths: " << m_nNumberOfSecondaryParticles << std::endl;
}


bool DetectorZPlaneMem::PermissionToResetPointer(int i_nStackSize)
{
	return true;
}
