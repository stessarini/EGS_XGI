# -*- coding: utf-8 -*-

###############################################################################
#
#   Python script - plots the MC simulated Talbot carpet
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
from import_detector_signal import get_simulation_time


print("############################################")

print("Import data...")

number_of_detectors=500
XPixels=30600
YPixels=3

all_signal=np.zeros([number_of_detectors,YPixels, XPixels])
for nSimCounter in range(number_of_detectors):
	binary_file_name='Carpet_G14d0p5_F41_20keV_Sim_' + str(nSimCounter+1) + '_detector1.bin'
	detector_signal=np.reshape(np.fromfile(binary_file_name,float),(YPixels,XPixels))
	all_signal[nSimCounter,:,:]=detector_signal[:,:]

del detector_signal




#get the carpet
middle_line_index=1
carpet=np.zeros([XPixels, number_of_detectors])
for nSimCounter in range(number_of_detectors):
	carpet[:,nSimCounter]=all_signal[nSimCounter, middle_line_index,:]

del all_signal




#crop carpet
FOV = 0.0102
window_size = 0.0022
window_size_in_pixels = int(XPixels * window_size / FOV)
cropped_carpet = carpet[int(XPixels / 2 - window_size_in_pixels / 2):int(XPixels / 2 + window_size_in_pixels / 2), :]
del carpet


#normalize each row of the carpet to have a max signal of 1
cropped_carpet_rownormalized = np.zeros([window_size_in_pixels, number_of_detectors])
for j in range(number_of_detectors):
	cropped_carpet_rownormalized[:,j] = cropped_carpet[:,j] / np.max(cropped_carpet[:,j])


fig_width = 11.0236 / 2
fig_height = 8.26772
my_dpi = 109


#rebin the MC SIGNAL
number_of_bins = int(cropped_carpet_rownormalized.shape[0]/10)
rebinned_carpet = np.zeros((number_of_bins, number_of_detectors))
for j in range(number_of_detectors):
	rebinned_carpet[:,j] = np.mean(cropped_carpet_rownormalized[:,j].reshape(number_of_bins, 10), 1)
	rebinned_carpet[:,j] = rebinned_carpet[:,j] / np.max(rebinned_carpet[:,j])



print("selected window of the carpet rebinned - normalization row wise - interpolation = None")
plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.imshow(rebinned_carpet * 255, cmap='gray', extent=[3.23,9.68,-11,11], aspect="auto", interpolation = 'None')
plt.xticks(np.around(np.linspace(3.23,9.68,5),decimals=2), fontsize = 25)
plt.yticks(np.around(np.linspace(-10,10,5),0),fontsize = 25)
plt.xlabel('Distance (cm)', fontsize = 27)
plt.ylabel(r'x ($\mu$m)', fontsize = 27,  labelpad=-28)
plt.subplots_adjust(bottom=0.13, top=0.97, right=0.92, left= 0.16)
plt.savefig("MCTalbotCarpet.pdf", dpi=300, format="pdf")
plt.show()



#simultion time

cpu_times=np.zeros(number_of_detectors)
elapsed_times=np.zeros(number_of_detectors)

for nSimCounter in range(number_of_detectors):
	log_file_name ='out_Carpet_G14d0p5_F41_20keV_Sim_' + str(nSimCounter+1) + '.txt'
	[t_cpu, t_ela] = get_simulation_time(log_file_name)
	cpu_times[nSimCounter] = t_cpu
	elapsed_times[nSimCounter] = t_ela

#plt.plot(cpu_times)
#plt.plot(elapsed_times)
#plt.show()

print("average cpu time: " + str(np.mean(cpu_times)))
print("average elapsed time: " + str(np.mean(elapsed_times)))
