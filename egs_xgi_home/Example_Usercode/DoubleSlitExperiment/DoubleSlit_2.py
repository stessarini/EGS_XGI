# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


from import_detector_signal import import_detector_signal_from_file

#XPixels=100000
#YPixels=3

FOV=0.0144

#Double slit parameters
b=2e-4;
a=2e-5;

#distance slits - detector
L=10.0;

#wavelength
lamb = (1.23984187541999*1e-10)/0.017;

#import signal double slit
file_name='DoubleSlit_d2em4_a2em5_CohPointSource_detector1'
[width, height, simulation_data] = import_detector_signal_from_file(file_name)
MC_signal_center_DS = simulation_data[int(round(height/2)), :]
del simulation_data

#x-axis
bin_width = FOV / width
x=np.linspace(-FOV/2 + bin_width/2.0, FOV/2 - bin_width / 2.0, width)

#plt.title("MC signal center")
#plt.plot(x + b / 2.0, MC_signal_center_DS)
#plt.grid()
#plt.show()


#import number of paths for each data point
file_name = 'DoubleSlit_d2em4_a2em5_CohPointSource_NOP_detector1.bin'
MC_number_of_paths = np.reshape(np.fromfile(file_name, dtype=np.uintc),(height, width))
MC_number_of_paths_center = MC_number_of_paths[int(round(height/2)), :]

#normalize MC signal
non_zero = np.greater(MC_number_of_paths_center,0)
MC_signal_center_DS_NOP2 = np.zeros(MC_number_of_paths_center.size)
MC_signal_center_DS_NOP2[non_zero] = MC_signal_center_DS[non_zero] / (MC_number_of_paths_center[non_zero]**2)


del MC_signal_center_DS, MC_number_of_paths_center, non_zero
total_MC_signal_center_DS_NOP2 = np.sum(MC_signal_center_DS_NOP2)







#the expected double slit signal
#cos_list_cor=np.cos(np.pi * (a+b) / lamb * (x + b/2.0)/np.sqrt(np.power(L,2) + np.power(x + b/2.0,2) ))
cos_list_cor=np.cos(np.pi * (a+b) / lamb * (x + b/2.0)/np.sqrt(L**2 + (x + b/2.0)**2))
#sinc_list_cor=np.pi * np.sinc( a / lamb * (x + b/2.0)/np.sqrt(np.power(L,2) + np.power(x + b/2.0,2) ))
sinc_list_cor=np.pi * np.sinc( a / lamb * (x + b/2.0)/np.sqrt(L**2 + (x + b/2.0)**2))

cos_2_list_cor=cos_list_cor**2
sinc_2_list_cor=sinc_list_cor**2

signal_thy_DS= cos_2_list_cor * sinc_2_list_cor
total_singal_thy_DS = np.sum(signal_thy_DS)






#normalize DS signals to same area with peak of thy signal = 1
print("normalized to same area")
double_slit_signal_thy = signal_thy_DS / total_singal_thy_DS
max__double_signal_thy = np.max(double_slit_signal_thy)
double_slit_signal_thy = double_slit_signal_thy / max__double_signal_thy
double_slit_signal_MC_NOP = MC_signal_center_DS_NOP2 / total_MC_signal_center_DS_NOP2 / max__double_signal_thy

del signal_thy_DS, MC_signal_center_DS_NOP2











#import blocked souble slit (BDS) results
file_name='DoubleSlit_d2em4_a2em5_CohPointSource_OneBlockedSlit_detector1'
[width, height, simulation_data] = import_detector_signal_from_file(file_name)
MC_signal_center_BDS = simulation_data[int(round(height/2)), :]
del simulation_data

file_name = 'DoubleSlit_d2em4_a2em5_CohPointSource_OneBlockedSlit_NOP_detector1.bin'
MC_number_of_paths_BDS = np.reshape(np.fromfile(file_name, dtype=np.uintc),(height, width))
MC_number_of_paths_center_BDS = MC_number_of_paths_BDS[int(round(height/2)), :]
del MC_number_of_paths_BDS




#normalize MC signal
non_zero = np.greater(MC_number_of_paths_center_BDS,0)
MC_signal_center_BDS_NOP2 = np.zeros(MC_number_of_paths_center_BDS.size)
MC_signal_center_BDS_NOP2[non_zero] = MC_signal_center_BDS[non_zero] / (MC_number_of_paths_center_BDS[non_zero]**2)


del MC_signal_center_BDS, MC_number_of_paths_center_BDS, non_zero
total_MC_signal_center_BDS_NOP2 = np.sum(MC_signal_center_BDS_NOP2)





#analythical single slit signal
deltaX = -(a +b) * 1.1 /2
sinc_list=np.pi * np.sinc( a / lamb * (x + b/2.0 - deltaX)/np.sqrt(np.power(L,2) + np.power(x + b/2.0 - deltaX,2) ))
signal_thy_SS = sinc_list**2

single_slit_signal_thy = signal_thy_SS / np.sum(signal_thy_SS)
max__single_slit_signal_thy = np.max(single_slit_signal_thy)
single_slit_signal_thy = single_slit_signal_thy / max__single_slit_signal_thy
blocked_double_slit_signal_MC_NOP = MC_signal_center_BDS_NOP2 / np.sum(MC_signal_center_BDS_NOP2) / max__single_slit_signal_thy

del signal_thy_SS, MC_signal_center_BDS_NOP2






#plot
print("normalized to area and peak height")
#plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
#plt.plot(10000 * (x + b/2.0 - deltaX), blocked_double_slit_signal_MC_NOP,'tab:orange', linestyle='solid', marker='', linewidth=2, label='Blocked slit')
#plt.plot(10000 * (x + b/2.0 - deltaX), single_slit_signal_thy, 'k', linestyle='dashed', linewidth=2)
#plt.plot(10000 * (x + b/2.0), double_slit_signal_MC_NOP,'m', linestyle='solid', marker='', linewidth=2, label='Double slit')
#plt.plot(10000 * (x + b/2.0), double_slit_signal_thy, 'k', linestyle='dashed', linewidth=2)
#plt.xlabel("x ($\mu$m)", fontsize = 27)
#plt.ylabel("Normalized detector signal", fontsize = 27)
#plt.xticks(np.arange(-60,80,20), fontsize = 25)
#plt.yticks(np.arange(0.0,1.1,0.20), fontsize = 25)
#plt.xlim((-70,70))
#plt.ylim((0.0,1.035))
#plt.grid()
##plt.legend(fontsize=25)
#plt.subplots_adjust(bottom=0.13, top=0.97, right=0.97, left= 0.11)
#plt.savefig("DoubleSlit_and_BlockedDoubleSlit.pdf", dpi=300, format="pdf")
#plt.show()

x_rebin4 = np.linspace(-FOV/2 + 4.0 * bin_width/2.0, FOV/2 - 4.0 * bin_width / 2.0, int(width/4))
x_rebin2 = np.linspace(-FOV/2 + 2.0 * bin_width/2.0, FOV/2 - 2.0 * bin_width / 2.0, int(width/2))
x_rebin8 = np.linspace(-FOV/2 + 8.0 * bin_width/2.0, FOV/2 - 8.0 * bin_width / 2.0, int(width/8))
#plt.plot(10000 * (x_rebin4 + b/2.0 - deltaX), np.mean(np.reshape(blocked_double_slit_signal_MC_NOP, (int(width/4),4)), axis=1))
#plt.plot(10000 * (x_rebin2 + b/2.0 - deltaX), np.mean(np.reshape(blocked_double_slit_signal_MC_NOP, (int(width/2),2)), axis=1))
#plt.plot(10000 * (x_rebin8 + b/2.0 - deltaX), np.mean(np.reshape(blocked_double_slit_signal_MC_NOP, (int(width/8),8)), axis=1))
#plt.plot(10000 * (x + b/2.0 - deltaX), single_slit_signal_thy)
#plt.show()

blocked_double_slit_signal_MC_NOP_rebinned = np.mean(np.reshape(blocked_double_slit_signal_MC_NOP, (int(width/8),8)), axis=1)

print("normalized to area and peak height")
#plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
#plt.plot(10000 * (x + b/2.0 - deltaX), blocked_double_slit_signal_MC_NOP_rebinned,'tab:orange', linestyle='solid', marker='', linewidth=2, label='Blocked slit')
#plt.plot(10000 * (x + b/2.0 - deltaX), single_slit_signal_thy, 'k', linestyle='dashed', linewidth=2)
#plt.plot(10000 * (x + b/2.0), double_slit_signal_MC_NOP,'m', linestyle='solid', marker='', linewidth=2, label='Double slit')
#plt.plot(10000 * (x + b/2.0), double_slit_signal_thy, 'k', linestyle='dashed', linewidth=2)
#plt.xlabel("x ($\mu$m)", fontsize = 27)
#plt.ylabel("Normalized detector signal", fontsize = 27)
#plt.xticks(np.arange(-60,80,20), fontsize = 25)
#plt.yticks(np.arange(0.0,1.1,0.20), fontsize = 25)
#plt.xlim((-70,70))
#plt.ylim((0.0,1.035))
#plt.grid()
##plt.legend(fontsize=25)
#plt.subplots_adjust(bottom=0.13, top=0.97, right=0.97, left= 0.11)
#plt.savefig("DoubleSlit_and_BlockedDoubleSlit.pdf", dpi=500, format="pdf")
#plt.show()


print("normalized to area and peak height")
plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
plt.plot(10000 * (x+ b/2.0 - deltaX), np.abs(single_slit_signal_thy - blocked_double_slit_signal_MC_NOP)/np.mean(blocked_double_slit_signal_MC_NOP),'tab:orange', linestyle='solid', marker='', linewidth=2, label='Blocked slit')
plt.plot(10000 * (x + b/2.0), np.abs(double_slit_signal_thy - double_slit_signal_MC_NOP)/np.mean(double_slit_signal_MC_NOP),'m', linestyle='solid', marker='', linewidth=2, label='Double slit')
plt.xlabel("x ($\mu$m)", fontsize = 27)
plt.ylabel(r'$\left|s_{thy} - s_{MC}\right| / \bar{s}_{thy}$', fontsize = 27)
#plt.ylabel(r'$\left|s_{thy} - s_{MC}\right| / \left<s_{thy}\right>$', fontsize = 27)
#plt.ylabel(r'$\epsilon / \mu$', fontsize = 27)
plt.xticks(np.arange(-60,80,20), fontsize = 25)
plt.yticks(np.arange(0.0,0.16,0.02), fontsize = 25)
#plt.yticks(np.arange(0.0,1.1,0.20), fontsize = 25)
plt.xlim((-70,70))
#plt.ylim((0.0,1.035))
plt.grid()
#plt.legend(fontsize=25)
plt.subplots_adjust(bottom=0.13, top=0.97, right=0.97, left= 0.11)
plt.savefig("DoubleSlit_and_BlockedDoubleSlit_differences.pdf", dpi=500, format="pdf")
plt.show()


single_slit_signal_thy_rebinned = np.mean(np.reshape(single_slit_signal_thy, (int(width/8),8)), axis=1)
double_slit_signal_MC_NOP_rebinned = np.mean(np.reshape(double_slit_signal_MC_NOP, (int(width/8),8)), axis=1)
double_slit_signal_thy_rebinned = np.mean(np.reshape(double_slit_signal_thy, (int(width/8),8)), axis=1)



print("normalized to area and peak height")
plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
plt.plot(10000 * (x_rebin8 + b/2.0 - deltaX), np.abs(single_slit_signal_thy_rebinned - blocked_double_slit_signal_MC_NOP_rebinned),'tab:orange', linestyle='solid', marker='', linewidth=2, label='Blocked slit')
plt.plot(10000 * (x_rebin8 + b/2.0), np.abs(double_slit_signal_thy_rebinned - double_slit_signal_MC_NOP_rebinned),'m', linestyle='solid', marker='', linewidth=2, label='Double slit')
plt.xlabel("x ($\mu$m)", fontsize = 27)
plt.ylabel("Absolute error", fontsize = 27)
plt.xticks(np.arange(-60,80,20), fontsize = 25)
#plt.yticks(np.arange(0.0,1.1,0.20), fontsize = 25)
plt.yticks(fontsize = 25)
plt.xlim((-70,70))
plt.ylim((0.0,0.011))
plt.grid()
#plt.legend(fontsize=25)
plt.subplots_adjust(bottom=0.13, top=0.97, right=0.97, left= 0.14)
plt.savefig("DoubleSlit_and_BlockedDoubleSlit_absolute_error_rebinned.pdf", dpi=500, format="pdf")
plt.show()




mean_single = np.mean(single_slit_signal_thy)
mean_double = np.mean(double_slit_signal_thy)
rmse_single = np.sqrt(np.mean((single_slit_signal_thy - blocked_double_slit_signal_MC_NOP)**2))
rmse_double = np.sqrt(np.mean((double_slit_signal_thy - double_slit_signal_MC_NOP)**2))

mean_single_rebinned = np.mean(single_slit_signal_thy_rebinned)
rmse_single_rebinned = np.sqrt(np.mean((single_slit_signal_thy_rebinned - blocked_double_slit_signal_MC_NOP_rebinned)**2))

print("rmse_single " + str(rmse_single))
print("mean_single/mean " + str(rmse_single / mean_single))
print("rmse_single_rebinned " + str(rmse_single_rebinned))
print("nrmse single rebinned " + str(rmse_single_rebinned / mean_single_rebinned))
print("rmse_double " + str(rmse_double))
print("rmse_double/mean " + str(rmse_double / mean_double))

print("rmse_single (%)" + str(100 * rmse_single))
print("rmse_double (%)" + str(100 * rmse_double))

print("normalized to area and peak height")
plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
plt.plot(10000 * (x_rebin8 + b/2.0 - deltaX), blocked_double_slit_signal_MC_NOP_rebinned,'tab:orange', linestyle='solid', marker='', linewidth=2, label='Blocked slit')
plt.plot(10000 * (x + b/2.0 - deltaX), single_slit_signal_thy, 'k', linestyle='dashed', linewidth=2)
plt.plot(10000 * (x + b/2.0), double_slit_signal_MC_NOP,'m', linestyle='solid', marker='', linewidth=2, label='Double slit')
plt.plot(10000 * (x + b/2.0), double_slit_signal_thy, 'k', linestyle='dashed', linewidth=2)
plt.xlabel("x ($\mu$m)", fontsize = 27)
plt.ylabel("Normalized detector signal", fontsize = 27)
plt.xticks(np.arange(-60,80,20), fontsize = 25)
plt.yticks(np.arange(0.0,1.1,0.20), fontsize = 25)
plt.xlim((-70,70))
plt.ylim((0.0,1.035))
plt.grid()
#plt.legend(fontsize=25)
plt.subplots_adjust(bottom=0.13, top=0.97, right=0.97, left= 0.11)
plt.savefig("DoubleSlit_and_BlockedDoubleSlit_reb.pdf", dpi=500, format="pdf")
plt.show()

non_zero_ss = np.greater(single_slit_signal_thy_rebinned,0)
err_rel_ss = np.zeros(single_slit_signal_thy_rebinned.size)
err_rel_ss = np.abs(single_slit_signal_thy_rebinned[non_zero_ss] - blocked_double_slit_signal_MC_NOP_rebinned[non_zero_ss]) / blocked_double_slit_signal_MC_NOP_rebinned[non_zero_ss]

non_zero_ds = np.greater(double_slit_signal_thy_rebinned,0)
err_rel_ds = np.zeros(double_slit_signal_thy_rebinned.size)
err_rel_ds = np.abs(double_slit_signal_thy_rebinned[non_zero_ds] - double_slit_signal_MC_NOP_rebinned[non_zero_ss]) / blocked_double_slit_signal_MC_NOP_rebinned[non_zero_ds]

print("rmse_single (%)" + str(100 * rmse_single))
print("rmse_double (%)" + str(100 * rmse_double))

print("normalized to area and peak height")
plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
plt.plot(10000 * (x_rebin8 + b/2.0 - deltaX), 100 * err_rel_ss,'tab:orange', linestyle='solid', marker='', linewidth=2, label='Blocked slit')
plt.plot(10000 * (x_rebin8 + b/2.0), 100 * err_rel_ds,'m', linestyle='solid', marker='', linewidth=2, label='Double slit')
plt.xlabel("x ($\mu$m)", fontsize = 27)
plt.ylabel("Relative error (%)", fontsize = 27)
plt.xticks(np.arange(-60,80,20), fontsize = 25)
#plt.yticks(np.arange(0.0,1.1,0.20), fontsize = 25)
plt.yticks(fontsize = 25)
plt.xlim((-70,70))
plt.ylim((0.0,10.035))
plt.grid()
#plt.legend(fontsize=25)
plt.subplots_adjust(bottom=0.13, top=0.97, right=0.97, left= 0.14)
plt.savefig("DoubleSlit_and_BlockedDoubleSlit_absolute_rel_error_rebinned.pdf", dpi=500, format="pdf")
plt.show()
