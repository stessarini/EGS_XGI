#include "xgi_global_variables.h"

EGS_Float phase_array[MXSTACK];
EGS_Float* e_fPhase = phase_array;

EGS_Float norm_array[MXSTACK];
EGS_Float* e_fLogNorm = norm_array;

bool primary_array[MXSTACK] = {false};
bool* e_bPrimary = primary_array;

const EGS_Float ec_fPi = 3.141592653589793238463;

const EGS_Float ec_fEnergyToWaveLength = 1.23984187541999*1e-10; /*\lambda[cm] = c_fEnergyToWaveLength/E[MeV]*/

const std::string ec_slibVersion = "EGS_XGI_v1.0";
