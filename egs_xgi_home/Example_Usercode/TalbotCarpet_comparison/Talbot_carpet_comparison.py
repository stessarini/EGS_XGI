###############################################################################
#
#   python script for analysis of Talbot carpet comparison
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
from import_detector_signal import get_simulation_time, import_detector_signal_from_file

Xbins = 30600
YBins = 3

FOVx = 0.0102
bin_width = FOVx/Xbins
x = np.linspace(-FOVx/2 + bin_width/2, FOVx/2 + bin_width/2, Xbins)

#crop window
window_size = 0.0022
window_size_in_pixels = int(Xbins * window_size / FOVx)

#rebin and average 10 consecutive bins
number_of_bins = int(window_size_in_pixels/10)

#smaple plane positions
index = np.arange(0,10) * 49


#Huygens simulated samples
N_G1_H = 2**np.arange(0,7) * 225
signals_Huygens = np.zeros((N_G1_H.size, index.size, number_of_bins))
t_cpu_H = np.zeros((N_G1_H.size, index.size))
for j in range(N_G1_H.size):
    for i in range(index.size):
        file_name = './Huygens_1e6_' + str(j) + '/TC_Hugens_1e6_exp_' + str(j) + "_ind_" + str(i) + '_detector1'
        [width, height, simulation_data] = import_detector_signal_from_file(file_name)
        cropped_signal = simulation_data[1,(int(Xbins/2 - window_size_in_pixels/2)):(int(Xbins/2 + window_size_in_pixels/2))]
        #plt.plot(cropped_signal)
        #plt.show()
        #print(cropped_signal.shape)
        rebinned_signal = np.mean(cropped_signal.reshape(number_of_bins,10),1)
        #normalize the singal w.r.t standard scalar product of N_G1_H.size * index.size dimensional vectors
        signals_Huygens[j,i,:] = rebinned_signal / (np.sum(rebinned_signal))
        log_file = './Huygens_1e6_' + str(j) + '/TC_Hugens_1e6_exp_' + str(j) + "_ind_" + str(i) + '.log'
        t_cpu_H[j,i] = get_simulation_time(log_file)[0]

signals_Huygens_5e5 = np.zeros((N_G1_H.size, index.size, number_of_bins))
t_cpu_H_5e5 = np.zeros((N_G1_H.size, index.size))
for j in range(N_G1_H.size):
    for i in range(index.size):
        file_name = './Huygens_5e5_' + str(j) + '/TC_Hugens_5e5_exp_' + str(j) + "_ind_" + str(i) + '_detector1'
        [width, height, simulation_data] = import_detector_signal_from_file(file_name)
        cropped_signal = simulation_data[1,(int(Xbins/2 - window_size_in_pixels/2)):(int(Xbins/2 + window_size_in_pixels/2))]
        #plt.plot(cropped_signal)
        #plt.show()
        #print(cropped_signal.shape)
        rebinned_signal = np.mean(cropped_signal.reshape(number_of_bins,10),1)
        #normalize the singal w.r.t standard scalar product of N_G1_H.size * index.size dimensional vectors
        signals_Huygens_5e5[j,i,:] = rebinned_signal / (np.sum(rebinned_signal))
        log_file = './Huygens_5e5_' + str(j) + '/TC_Hugens_5e5_exp_' + str(j) + "_ind_" + str(i) + '.log'
        t_cpu_H_5e5[j,i] = get_simulation_time(log_file)[0]


#Fourier splitting simulated samples
N_G1_F = np.arange(5,43,2)
signals_Fourier_1e6 = np.zeros((N_G1_F.size, index.size, number_of_bins))
n = 0
t_cpu_F_1e6 = np.zeros((N_G1_F.size, index.size))
for j in N_G1_F:
    for i in range(index.size):
        file_name = './Fourier_1e6_' + str(j) + '/TC_Fourier_1e6_NG1_' + str(j) + "_ind_" + str(i) + '_detector1'
        [width, height, simulation_data] = import_detector_signal_from_file(file_name)
        cropped_signal = simulation_data[1,(int(Xbins/2 - window_size_in_pixels/2)):(int(Xbins/2 + window_size_in_pixels/2))]
        rebinned_signal = np.mean(cropped_signal.reshape(number_of_bins,10),1)
        #plt.plot(rebinned_signal)
        #plt.show()
        signals_Fourier_1e6[n,i,:] = rebinned_signal / (np.sum(rebinned_signal))
        log_file = './Fourier_1e6_' + str(j) + '/TC_Fourier_1e6_NG1_' + str(j) + "_ind_" + str(i) + '.log'
        t_cpu_F_1e6[n,i] = get_simulation_time(log_file)[0]
    n += 1


#Fourier splitting simulated samples
signals_Fourier_4e6 = np.zeros((N_G1_F.size, index.size, number_of_bins))
n = 0
t_cpu_F_4e6 = np.zeros((N_G1_F.size, index.size))
for j in N_G1_F:
    for i in range(index.size):
        file_name = './Fourier_4e6_' + str(j) + '/TC_Fourier_4e6_NG1_' + str(j) + "_ind_" + str(i) + '_detector1'
        [width, height, simulation_data] = import_detector_signal_from_file(file_name)
        cropped_signal = simulation_data[1,(int(Xbins/2 - window_size_in_pixels/2)):(int(Xbins/2 + window_size_in_pixels/2))]
        rebinned_signal = np.mean(cropped_signal.reshape(number_of_bins,10),1)
        #plt.plot(rebinned_signal)
        #plt.show()
        signals_Fourier_4e6[n,i,:] = rebinned_signal / (np.sum(rebinned_signal))
        log_file = './Fourier_4e6_' + str(j) + '/TC_Fourier_4e6_NG1_' + str(j) + "_ind_" + str(i) + '.log'
        t_cpu_F_4e6[n,i] = get_simulation_time(log_file)[0]
    n += 1

signals_Fourier_8e6 = np.zeros((N_G1_F.size, index.size, number_of_bins))
n = 0
t_cpu_F_8e6 = np.zeros((N_G1_F.size, index.size))
for j in N_G1_F:
    for i in range(index.size):
        file_name = './Fourier_8e6_' + str(j) + '/TC_Fourier_8e6_NG1_' + str(j) + "_ind_" + str(i) + '_detector1'
        [width, height, simulation_data] = import_detector_signal_from_file(file_name)
        cropped_signal = simulation_data[1,(int(Xbins/2 - window_size_in_pixels/2)):(int(Xbins/2 + window_size_in_pixels/2))]
        rebinned_signal = np.mean(cropped_signal.reshape(number_of_bins,10),1)
        #plt.plot(rebinned_signal)
        #plt.show()
        signals_Fourier_8e6[n,i,:] = rebinned_signal / (np.sum(rebinned_signal))
        log_file = './Fourier_8e6_' + str(j) + '/TC_Fourier_8e6_NG1_' + str(j) + "_ind_" + str(i) + '.log'
        t_cpu_F_8e6[n,i] = get_simulation_time(log_file)[0]
    n += 1



#Wave propagator samples
all_signal_wp = np.load("../WP_patterns/WP_rebinned_carpet.npy")
signal_WP = np.zeros((index.size, number_of_bins))
signal_WP = np.transpose(all_signal_wp[:,index])
for j in range(index.size):
    signal_WP[j,:] = signal_WP[j,:] / (np.sum(signal_WP[j,:]))

#calculate difference measures


mean_WP = np.mean(signal_WP)
mean_H = np.mean(signals_Huygens, (1,2))
mean_F_4e6 = np.mean(signals_Fourier_4e6, (1,2))

var_WP = np.var(signal_WP)
var_H = np.var(signals_Huygens, (1,2))
var_F_1e6 = np.var(signals_Fourier_1e6, (1,2))
var_F_4e6 = np.var(signals_Fourier_4e6, (1,2))
var_F_8e6 = np.var(signals_Fourier_8e6, (1,2))




#calculate Pearson correlations:
cor_H = np.zeros(N_G1_H.size)
for j in range(N_G1_H.size):
    cov = np.cov(signal_WP.flatten(),signals_Huygens[j,:,:].flatten())
    #print(cov)
    cor_H[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
#print("cor_H: " + str(cor_H))

cor_H_5e5 = np.zeros(N_G1_H.size)
for j in range(N_G1_H.size):
    cov = np.cov(signal_WP.flatten(),signals_Huygens_5e5[j,:,:].flatten())
    #print(cov)
    cor_H_5e5[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
#print("cor_H_5e5: " + str(cor_H))


cor_F_8e6 = np.zeros(N_G1_F.size)
for j in range(N_G1_F.size):
        cov = np.cov(signal_WP.flatten(), signals_Fourier_8e6[j,:,:].flatten())
        cor_F_8e6[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
#rint("cor_F_8e6: " + str(cor_F_8e6))


cor_F_4e6 = np.zeros(N_G1_F.size)
for j in range(N_G1_F.size):
        cov = np.cov(signal_WP.flatten(), signals_Fourier_4e6[j,:,:].flatten())
        cor_F_4e6[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
#print("cor_F_4e6: " + str(cor_F_4e6))

cor_F_1e6 = np.zeros(N_G1_F.size)
for j in range(N_G1_F.size):
        cov = np.cov(signal_WP.flatten(), signals_Fourier_1e6[j,:,:].flatten())
        cor_F_1e6[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
#print("cor_F_1e6: " + str(cor_F_1e6))


mean_cpu_time_H = np.mean(t_cpu_H,1)
mean_cpu_time_H_5e5 = np.mean(t_cpu_H_5e5,1)
mean_cpu_time_F_1e6 = np.mean(t_cpu_F_1e6,1)
mean_cpu_time_F_4e6 = np.mean(t_cpu_F_4e6,1)
mean_cpu_time_F_8e6 = np.mean(t_cpu_F_8e6,1)

print("relative cpu time for last (most accurate) cases: " + str(mean_cpu_time_H[-1] / mean_cpu_time_F_4e6[-1]))
print("correlation coeff:")
print(cor_H[-1])
print(cor_F_4e6[-1])
print(" relative correlation coeff")
print(cor_H[-1] / cor_F_4e6[-1])



print("++++++++")
print("get relative cpu time closest case to best Huygens case: ")
d_corr_1e6 = (np.abs(cor_H[-1] / cor_F_1e6[:] - 1))
d_corr_4e6 = (np.abs(cor_H[-1] / cor_F_4e6[:] - 1))
d_corr_8e6 = (np.abs(cor_H[-1] / cor_F_8e6[:] - 1))

argmin_1e6 = np.argmin(d_corr_1e6)
argmin_4e6 = np.argmin(d_corr_4e6)
argmin_8e6 = np.argmin(d_corr_8e6)

print("mean cpu time best Huygens splitting case:")
print(mean_cpu_time_H[-1])
print("cor best Huygens splitting case:")
print(cor_H[-1])
print("++++++++")
print("F_1e6:")
print("argmin:")
print(argmin_1e6)
print("differece score:")
print(d_corr_1e6[argmin_1e6])
print("cor F:")
print(cor_F_1e6[argmin_1e6])
print("NG1")
print(N_G1_F[argmin_1e6])
print("time imporvement:")
print(mean_cpu_time_H[-1] / mean_cpu_time_F_1e6[argmin_1e6])
print("average simulation time")
print(mean_cpu_time_F_1e6[argmin_1e6])

print("++++++++")
print("F_4e6:")
print("argmin:")
print(argmin_4e6)
print("differece score:")
print(d_corr_4e6[argmin_4e6])
print("cor F:")
print(cor_F_4e6[argmin_4e6])
print("NG1")
print(N_G1_F[argmin_4e6])
print("time imporvement:")
print(mean_cpu_time_H[-1] / mean_cpu_time_F_4e6[argmin_4e6])
print("average simulation time")
print(mean_cpu_time_F_4e6[argmin_4e6])

print("++++++++")
print("F_8e6:")
print("argmin:")
print(argmin_8e6)
print("differece score:")
print(d_corr_8e6[argmin_8e6])
print("cor F:")
print(cor_F_8e6[argmin_8e6])
print("NG1")
print(N_G1_F[argmin_8e6])
print("time imporvement:")
print(mean_cpu_time_H[-1] / mean_cpu_time_F_8e6[argmin_8e6])
print("average simulation time")
print(mean_cpu_time_F_8e6[argmin_8e6])
print("++++++++")








#plot all correlation coefficients

plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
plt.plot(mean_cpu_time_H_5e5, cor_H_5e5, label=r'H: $N_H= 5\times10^5$', linewidth=2, color='k')
plt.plot(mean_cpu_time_H, cor_H, label=r'H: $N_H=10^6$', linewidth=2, color='m')
plt.plot(mean_cpu_time_F_1e6, cor_F_1e6, label=r'F: $N_H=10^6$', linewidth=2, color='green')
plt.plot(mean_cpu_time_F_4e6, cor_F_4e6, label=r'F: $N_H=4\times10^6$', linewidth=2, color='tab:orange')
plt.plot(mean_cpu_time_F_8e6, cor_F_8e6, label=r'F: $N_H=8\times10^6$', linewidth=2, color='blue')
plt.xscale('log')
plt.xlabel(r'$t_{cpu}$ (s)' , fontsize = 27)
plt.ylabel(r'Correlation coefficient', fontsize = 27)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
#plt.xlim((-70,70))
#plt.ylim((-70,70))
plt.legend(fontsize = 25, loc = 'lower right')
plt.grid()
#plt.subplots_adjust(bottom=0.12, top=0.97, right=0.98, left= 0.11)
plt.subplots_adjust(bottom=0.13, top=0.97, right=0.98, left= 0.12)
plt.savefig("cor_tcpu_Huygens_Fourier_10carpetPositions_NH_ext_3_lab.pdf", dpi=500, format="pdf")
plt.show()

mean_cpu_time_F = np.array([mean_cpu_time_F_1e6, mean_cpu_time_F_4e6])
print(mean_cpu_time_F.size)
cor_F = np.array([cor_F_1e6, cor_F_4e6])
