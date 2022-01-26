# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from import_detector_signal import get_simulation_time


print("############################################")
print("python for creating Talbot carpet from simulaiton output")
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


plt.plot(cropped_carpet_rownormalized[:,0])
plt.show()

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
#plt.ylabel((u'x (\u03BCm)'), fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
plt.ylabel(r'x ($\mu$m)', fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
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

plt.plot(cpu_times)
plt.plot(elapsed_times)
plt.show()

print("average cpu time: " + str(np.mean(cpu_times)))
print("average elapsed time: " + str(np.mean(elapsed_times)))





#All that follows can be removed





number_of_paths=np.zeros([number_of_detectors,YPixels, XPixels])
for nSimCounter in range(number_of_detectors):
	binary_file_name='Carpet_G14d0p5_F41_20keV_Sim_' + str(nSimCounter+1) + '_NOP_detector1.bin'
	detector_signal=np.reshape(np.fromfile(binary_file_name,dtype=np.uintc),(YPixels,XPixels))
	number_of_paths[nSimCounter,:,:]=detector_signal[:,:]


del detector_signal

NP_carpet=np.zeros([XPixels, number_of_detectors])
for nSimCounter in range(number_of_detectors):
	NP_carpet[:,nSimCounter]=number_of_paths[nSimCounter, middle_line_index,:]

del number_of_paths

cropped_NP_carpet = NP_carpet[int(XPixels / 2 - window_size_in_pixels / 2):int(XPixels / 2 + window_size_in_pixels / 2), :]
del NP_carpet

#print("selected window of the carpet rebinned - normalization row wise - interpolation = None")
#plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
#plt.imshow(cropped_NP_carpet/np.max(cropped_NP_carpet) * 255, cmap='gray', extent=[2.74,8.23,-11,11], aspect="auto", interpolation = 'None')
#plt.xticks(np.around(np.linspace(2.74,8.23,5),decimals=2), fontsize = 25)
#plt.yticks(np.around(np.linspace(-10,10,5),0),fontsize = 25)
#plt.xlabel('Distance (cm)', fontsize = 27)
##plt.ylabel((u'x (\u03BCm)'), fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
#plt.ylabel("x ($\mu$m)", fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
#plt.subplots_adjust(bottom=0.13, top=0.97, right=0.92, left= 0.16)
#plt.show()

non_zero = np.greater(cropped_NP_carpet,0)
carpet_NOP = np.zeros(non_zero.shape)
carpet_NOP[non_zero] =  cropped_carpet[non_zero] / (cropped_NP_carpet[non_zero]**2)
for j in range(number_of_detectors):
	carpet_NOP[:,j] = carpet_NOP[:,j] / np.max(carpet_NOP[:,j])

rebinned_carpet_NOP = np.zeros((int(number_of_bins), number_of_detectors))
for j in range(number_of_detectors):
	rebinned_carpet_NOP[:,j] = np.mean(carpet_NOP[:,j].reshape(int(number_of_bins), 10), 1)
	rebinned_carpet_NOP[:,j] = rebinned_carpet_NOP[:,j] / np.max(rebinned_carpet_NOP[:,j])

plt.plot(rebinned_carpet_NOP[:,0])
plt.show()

plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.imshow(rebinned_carpet_NOP * 255, cmap='gray', extent=[3.23,9.68,-11,11], aspect="auto", interpolation = 'None')
plt.xticks(np.around(np.linspace(3.23,9.68,5),decimals=2), fontsize = 25)
plt.yticks(np.around(np.linspace(-10,10,5),0),fontsize = 25)
plt.xlabel('Distance (cm)', fontsize = 27)
#plt.ylabel((u'x (\u03BCm)'), fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
plt.ylabel("x ($\mu$m)", fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
plt.subplots_adjust(bottom=0.13, top=0.97, right=0.92, left= 0.16)
plt.savefig("MCTalbotCarpet_NOP.pdf", dpi=300, format="pdf")
plt.show()




all_signal_wp = np.load("/media/gerhard/SeagateExpansionDrive/PhD/MC/EGSnrc18/EGS_XGI_pub2/egs_xgi_home/Parameter_tests/TalbotCarpets/WP_TC_python/WP_rebinned_carpet.npy")

#set the signals to same average
print("mean:")
print(np.mean(rebinned_carpet_NOP))
print(np.mean(all_signal_wp))
carpet_renormalized = rebinned_carpet_NOP / np.mean(rebinned_carpet_NOP) * np.mean(all_signal_wp)
print(np.mean(carpet_renormalized))


plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.imshow(np.abs(carpet_renormalized - all_signal_wp) * 255, cmap='gray', vmin=0, vmax=255, extent=[3.23,9.68,-11,11], aspect="auto", interpolation = 'None')
plt.xticks(np.around(np.linspace(3.23,9.68,5),decimals=2), fontsize = 25)
plt.yticks(np.around(np.linspace(-10,10,5),0),fontsize = 25)
plt.xlabel('Distance (cm)', fontsize = 27)
#plt.ylabel((u'x (\u03BCm)'), fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
plt.ylabel("x ($\mu$m)", fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
plt.subplots_adjust(bottom=0.13, top=0.97, right=0.92, left= 0.16)
plt.savefig("MCTalbotCarpet_NOP_diff_WP.pdf", dpi=300, format="pdf")
plt.show()



plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.imshow(carpet_renormalized * 255, cmap='gray', extent=[3.23,9.68,-11,11],vmin=0, vmax=255, aspect="auto", interpolation = 'None')
plt.imshow(np.abs(carpet_renormalized - all_signal_wp) * 255, cmap='Purples',alpha=0.5, extent=[3.23,9.68,-11,11], aspect="auto", interpolation = 'None')
plt.xticks(np.around(np.linspace(3.23,9.68,5),decimals=2), fontsize = 25)
plt.yticks(np.around(np.linspace(-10,10,5),0),fontsize = 25)
plt.xlabel('Distance (cm)', fontsize = 27)
#plt.ylabel((u'x (\u03BCm)'), fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
plt.ylabel("x ($\mu$m)", fontsize = 27,  labelpad=-28)#u'x [\u03BCm]'
plt.subplots_adjust(bottom=0.13, top=0.97, right=0.92, left= 0.16)
#plt.savefig("MCTalbotCarpet_NOP_diff_WP.pdf", dpi=300, format="pdf")
plt.show()

print("maxima")
print(np.max(carpet_renormalized))
print(np.max(all_signal_wp))
rmse = np.sqrt(np.mean((carpet_renormalized - all_signal_wp)**2))
nrmse = rmse / np.mean(all_signal_wp)
print("rmse = " + str(rmse))
print("nrmse = " + str(nrmse))
print("rmse in grey values= " + str(255 * rmse))
[N1,N2] = carpet_renormalized.shape
print([N1,N2])
print("rmse alt: " + str(np.sqrt(np.sum((carpet_renormalized - all_signal_wp)**2)/N1/N2)))


max_diff = (np.max((carpet_renormalized - all_signal_wp)**2))
print("max_diff: " + str(max_diff))

var = np.var(carpet_renormalized.flatten() - all_signal_wp.flatten())
print("var: " + str(var))

SNR = np.mean(carpet_renormalized) / np.sqrt(var)
print("SNR = " + str(SNR))

print(1/rmse)

cov = np.cov(carpet_renormalized.flatten(), all_signal_wp.flatten())
cor = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
print("cor = " + str(cor))
print(np.corrcoef(carpet_renormalized.flatten(), all_signal_wp.flatten()))
