# -*- coding: utf-8 -*-
###############################################################################
#
#   python analysis file for lab setup simulation
#   Copyright (C) 2020  ETH Zürich
#
#   This file is part of the EGS_XGI - an X-ray grating interferometry
#   extension for EGSnrc.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as published
#   by the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
###############################################################################
#
#   Author:     Stefan Tessarini
#
#
#
###############################################################################
import numpy as np
import matplotlib.pyplot as plt


from import_detector_signal import import_detector_signal_from_file, calculate_phase_stepping_curve, get_visibility_map



pixelWidth = 0.0075
FOV = 0.5025

nPixelX = int(np.floor(FOV / pixelWidth))
print("number of pixels in x-direction: " + str(nPixelX))
nPixelY = 1;


print("############################################")
print("Import data...")
nNumberOfSimulations = 10

binary_file_name='LabSetup0_detectorMem1'
[width, height, simulation_data] = import_detector_signal_from_file(binary_file_name)

for j in range(1,10):
    binary_file_name='LabSetup' + str(j) + '_detectorMem1'
    print("load: " + binary_file_name)
    [width, height, simulation_data_2] = import_detector_signal_from_file(binary_file_name)
    simulation_data = simulation_data + simulation_data_2

del simulation_data_2



#####################################
#imprt reference without graings
binary_file_name='LabSetup_NoG_0_detectorMem1'
[width, height, simulation_data_NoG] = import_detector_signal_from_file(binary_file_name)
for j in range(1,10):
    binary_file_name='LabSetup_NoG_' + str(j) + '_detectorMem1'
    print("load: " + binary_file_name)
    [width, height, simulation_data_2] = import_detector_signal_from_file(binary_file_name)
    simulation_data_NoG += simulation_data_2


del simulation_data_2
averaged_signal_NoG =sum(simulation_data_NoG)

######################################



#rebin no grating signal into 67 pixel signals (neglect the 2 boundary pixels)
M = 40000
rebinned_average_MC_intensity_NOG = np.sum(np.reshape(averaged_signal_NoG, (int(width/(M)), (M))), axis = 1)[1:65]




normalization_constant = np.mean(rebinned_average_MC_intensity_NOG)


######################################


total_intensity_NG0 = np.sum(simulation_data_NoG)
total_intensity = np.sum(simulation_data)
print("--------------------------------------------------------------------------")
print("relative MC signal: ('MC signal'/'MC signal no gratings') " + str(total_intensity / total_intensity_NG0))

###############################################################
#spectrum:
energies = np.array([0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041])
probabilities = np.array([4.437995661681658904e-02,4.764742814870436688e-02,4.951928229267677645e-02,5.039431363363445354e-02,5.052983537739734465e-02,5.012325928304973977e-02,4.932991172646215838e-02,4.824028563960389149e-02,4.693288908610587090e-02,4.545907408060136101e-02,4.385757275001288952e-02,4.216006168266742127e-02,4.038555493014862885e-02,3.855087942464249795e-02,3.666867102173651527e-02,3.474917771964486368e-02,3.279963008751592135e-02,3.082539474740627478e-02,2.883111256904013694e-02,2.682049261697410855e-02,2.479643213260912873e-02,2.276149851103172542e-02,2.071727722716266984e-02,1.866533447707846224e-02,1.660692428997067108e-02,1.454310286062679772e-02,1.247471706611745848e-02,1.040243044455944552e-02,8.326830858456129442e-03,6.248410015625116033e-03,4.167576966187052333e-03,2.084681715733581030e-03])


#attenuation coefficients
mu_SI = np.array([77.489250971352945, 58.434056900334653, 45.179366702956095, 35.595827502772558, 28.539438401833078, 23.246034894745648, 19.160729788534763, 15.989773088615539, 13.479732420725533, 11.473387262596797, 9.8542424227788743, 8.5241336730406356, 7.4300881022300738, 6.5219598179508091, 5.7623080090096837, 5.1191859083167719, 4.5713422683235709, 4.1039372510164993, 3.7028887809914361, 3.3567300988167235, 3.0565669576836458, 2.7929737449358671, 2.5624974581941502, 2.3599998547099306, 2.1814279879276222, 2.0230512220274486, 1.881946889667339, 1.7562374718610934, 1.6438961080448327, 1.5431686291222797, 1.4526093651159406, 1.3705389085813759])
mu_AIR = np.array([0.0059021352785648554, 0.0044287402370836474, 0.0034178032797347688, 0.0026989039826863697, 0.0021762664677626542, 0.0017883519077147112, 0.0014932017149677369, 0.0012662494491440064, 0.0010883980431343437, 0.00094747763015618836, 0.00083461802611759059, 0.00074266204211392037, 0.0006675150299284944, 0.00060551915990923375, 0.00055394624872729, 0.0005105611166340069, 0.00047379686507023204, 0.00044256635975329717, 0.00041585808281287588, 0.0003928620446498325, 0.00037298299948653959, 0.00035557234884463082, 0.00034037042659868498, 0.00032701746065019418, 0.00031524332242772595, 0.00030481088978573162, 0.00029552288372937613, 0.00028722826484503983, 0.00027979426481708909, 0.0002730977333741131, 0.00026706277649997607, 0.00026156124188115811])
mu_AU = np.array([2188.658911862, 1707.8184842863566, 3469.9017962917937, 2795.8413364694693, 3189.4935892506769, 3101.3558926876244, 2628.9725146786459, 2251.0031615937969, 1942.536588649883, 1689.0047369949643, 1479.120054073167, 1301.8522818990036, 1152.6818101014082, 1025.2305509056896, 916.34581925247312, 822.34311378232462, 740.75033183791788, 669.9029537725753, 608.07993597890788, 553.84595716708634, 506.0703899953441, 463.45400634336499, 425.63918463726935, 391.94047288921797, 361.81490845088121, 334.73577268981546, 310.29137362387524, 288.24124322292312, 268.29680552336094, 250.20609210648257, 233.75212658722012, 218.66380826321745])
mu_IR = np.array([2339.67, 1824.85, 3728.45, 4171.03, 4006.42, 3362.7, 2847.19, 2435.19, 2098.34, 1821.49, 1592.65, 1399.88, 1237.81, 1100.54, 983.437, 882.325, 794.553, 718.375, 651.906, 593.634, 542.278, 496.505, 455.874, 419.689, 387.338, 358.27, 332.022, 308.349, 286.938, 267.526, 249.882, 233.692])

###############################################################
#Calculate classically (Beer-Lambert) expected total signal reduction from gratings

#transmission through absorbing G0 given by duty cycle
fDutyCycle_G0 = 0.5;

#transmission through G1
fDutyCycle = 0.5
periodG1 = 0.00015
thickness_G1 = 0.0025

NormTransmissionFunctionA = np.exp(-thickness_G1 * mu_SI / 2.0)
NormTransmissionFunctionB = np.exp(-thickness_G1 * mu_AIR / 2.0)

average_G1_attenuation = fDutyCycle * NormTransmissionFunctionA**2 + (1.0 - fDutyCycle) * NormTransmissionFunctionB**2


#transmission through air and silicon
d_SI_tot = 0.0220 + 0.0225
d_AIR_tot = 90.4 - d_SI_tot
abs_SI = np.exp(-mu_SI * d_SI_tot)
abs_AIR = np.exp(-mu_AIR * d_AIR_tot)

#An estimation for the signals (over an area significantly larger than the grating periods) without and with gratings:
expected_total_signal_NOG = np.sum(probabilities * abs_SI * abs_AIR)

expected_total_signal_with_gratings = fDutyCycle_G0 * np.sum(probabilities * abs_SI * abs_AIR * average_G1_attenuation)

classically_expected_intensity_reduction = expected_total_signal_with_gratings / expected_total_signal_NOG
print("classically_expected_intensity_reduction" + str(classically_expected_intensity_reduction))




print("MC over expected rel. signals: " + str(total_intensity / total_intensity_NG0 / classically_expected_intensity_reduction))
print("Difference in relative signals: " + str(1 - total_intensity / total_intensity_NG0 / classically_expected_intensity_reduction))

print("--------------------------------------------------------------------------")

##############################################################
#plot the MC signals with and without gratings
averaged_signal=sum(simulation_data)
rebinned_average_MC_intensity_G0G1 = np.sum(np.reshape(averaged_signal, (int(width/10), 10)), axis = 1)
x_10=np.linspace(-FOV/2.0,FOV/2.0,int(width/10))
x_M = np.linspace(-FOV/2.0,FOV/2.0,int(width/(M)) )


##########################3
# Plot parameters
#
fig_width = 11.0236
fig_height = 8.26772
my_dpi = 109


#pixel centers:
x_pixel = np.linspace(-FOV / 2.0 + pixelWidth / 2.0, FOV / 2.0 - pixelWidth / 2.0, nPixelX)
middle = int(np.floor(nPixelX / 2))
x_middle = x_pixel[middle]

pixel_borderX_0 = 10000 * np.array([x_pixel[middle - 1] - pixelWidth / 2.0, x_pixel[middle - 1] - pixelWidth / 2.0])
pixel_borderX_1 = 10000 * np.array([x_middle - pixelWidth / 2.0, x_middle - pixelWidth / 2.0])
pixel_borderX_2 = 10000 * np.array([x_pixel[middle + 1] - pixelWidth / 2.0, x_pixel[middle + 1] - pixelWidth / 2.0])
pixel_borderX_3 = 10000 * np.array([x_pixel[middle + 2] - pixelWidth / 2.0, x_pixel[middle + 2] - pixelWidth / 2.0])
pixel_borderY = np.array([0.0, 1.1])
plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.plot(10000 * x_M[1:65], rebinned_average_MC_intensity_NOG / normalization_constant, color='r', marker="", linestyle='solid', linewidth=2)
plt.plot(10000 * x_10, rebinned_average_MC_intensity_G0G1 / normalization_constant * M/10, color='b', marker="", linestyle='solid', linewidth=2)
plt.plot(pixel_borderX_0, pixel_borderY, color='grey', marker="", linestyle="--", linewidth = 2)
plt.plot(pixel_borderX_1, pixel_borderY, color='grey', marker="", linestyle="--", linewidth = 2)
plt.plot(pixel_borderX_2, pixel_borderY, color='grey', marker="", linestyle="--", linewidth = 2)
plt.plot(pixel_borderX_3, pixel_borderY, color='grey', marker="", linestyle="--", linewidth = 2)
plt.xticks(np.array([-112.5, -75, -37.5, 0.0, 37.5, 75, 112.5]), fontsize=25)
plt.yticks(fontsize=25)
plt.ylim(0,1.1)
plt.xlim((10000 * (x_middle - 3 * pixelWidth / 2.0) - 10, 10000 * (x_middle + 3.0 * pixelWidth / 2.0) + 10))
plt.xlabel('x($\mu$m)', fontsize=27)
plt.ylabel('Normalized MC signal', fontsize=27)
plt.grid()
plt.subplots_adjust(bottom=0.13, top=0.957, right=0.95, left= 0.155)
plt.savefig("LabSetup_results.pdf", dpi=300, format="pdf")
plt.show()


print("--------------------------------------------------------------------------")


###########################3
#get the visibility from phase stepping curve
nPhaseSteps = 10 #should divide 2*nMask
periodG2 = 0.0003
nMask = int(np.round(width * periodG2 / FOV /2 ))
print("nMask = " + str(nMask))
number_of_periods = int(np.floor(width / (2 * nMask)))
datapoints_per_pixel = int(np.floor(width * pixelWidth / FOV))
cropped_width = datapoints_per_pixel * nPixelX;
to_crop = width - cropped_width
print("to_crop: " + str(to_crop))
cropped_simulation_data = simulation_data[:, int(np.floor(to_crop/2)):(width-int(np.floor(to_crop/2)))]
print("size of cropped data: " + str(cropped_simulation_data.shape))

simulated_psc = calculate_phase_stepping_curve(nPixelX, nPixelY, nMask, nPhaseSteps, cropped_simulation_data)


print("-------------------------------------------------------")

simulated_visibility_map_2NG0 = get_visibility_map(simulated_psc)
print("-------------------------------------------------------")

print("simulated visibility in the center: " + str(simulated_visibility_map_2NG0[0,int(nPixelX/2)]))
print("average visibility: " + str(np.mean(simulated_visibility_map_2NG0[0,:])))
print("average visibility without first and last pixel: " + str(np.mean(simulated_visibility_map_2NG0[0,1:65])))
print("maximum visibility: " + str(np.max(simulated_visibility_map_2NG0[0,:])))
print("maximum visibility without first and last pixel: " + str(np.max(simulated_visibility_map_2NG0[0,1:65])))
print("minimum visibility: " + str(np.min(simulated_visibility_map_2NG0[0,:])))
print("minimum visibility without first and last pixel: " + str(np.min(simulated_visibility_map_2NG0[0,1:65])))




######################################
#visibility correction
#The visibility above is caluclated with a binary grating
#The leakage thrugh the absorber and the attenuation of the filling can reduce Visibility

#for a monochromatic beam the correction factor amounts to
print("-------------------------------------------------------")

print("visibility correction G2 leakage for monochromatic beams. Consider mean energy")
ID_char = 11
thickness_G2 = 0.003
print("mean energy " + str(energies[ID_char]))
print("transmission AU (from EGSnrc): " + str(np.exp(-mu_AU[ID_char] * thickness_G2)))
print("transmission SI: " + str(np.exp(-mu_SI[ID_char] * thickness_G2)))

v_corr_monpo_AUSI_char = (np.exp(-mu_SI[ID_char] * thickness_G2) - np.exp(-mu_AU[ID_char] * thickness_G2)) / (np.exp(-mu_SI[ID_char] * thickness_G2) + np.exp(-mu_AU[ID_char] * thickness_G2))

print("visibility correction mono AU-SI grating: " + str(v_corr_monpo_AUSI_char))
print("visibility corrected, mono AU-SI grating: " + str(v_corr_monpo_AUSI_char * np.mean(simulated_visibility_map_2NG0[0,1:65])))



print("-------------------------------------------------------")
print("averaged correction factor over all energies")
av_corr = np.sum(probabilities * (np.exp(-mu_SI * thickness_G2) - np.exp(-mu_AU * thickness_G2)) / (np.exp(-mu_SI * thickness_G2) + np.exp(-mu_AU * thickness_G2)))
print("visibility correction poly AU-SI grating: " + str(av_corr))
print("visibility corrected poly AU-SI grating: " + str(av_corr * np.mean(simulated_visibility_map_2NG0[0,1:65])))
