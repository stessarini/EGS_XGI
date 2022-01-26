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
print("N_G1_H = " + str(N_G1_H))
for j in range(N_G1_H.size):
    for i in range(index.size):
        file_name = './Huygens_' + str(j) + '/TC_Hugens_exp_' + str(j) + "_ind_" + str(i) + '_detector1'
        [width, height, simulation_data] = import_detector_signal_from_file(file_name)
        cropped_signal = simulation_data[1,(int(Xbins/2 - window_size_in_pixels/2)):(int(Xbins/2 + window_size_in_pixels/2))]
        #plt.plot(cropped_signal)
        #plt.show()
        #print(cropped_signal.shape)
        rebinned_signal = np.mean(cropped_signal.reshape(number_of_bins,10),1)
        #normalize the singal w.r.t standard scalar product of N_G1_H.size * index.size dimensional vectors
        signals_Huygens[j,i,:] = rebinned_signal / (np.sum(rebinned_signal))
        log_file = './Huygens_' + str(j) + '/TC_Hugens_exp_' + str(j) + "_ind_" + str(i) + '.log'
        t_cpu_H[j,i] = get_simulation_time(log_file)[0]
    print(np.sum(signals_Huygens[j,:,:]))

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
    print(np.sum(signals_Huygens[j,:,:]))


#Fourier splitting simulated samples
N_G1_F = np.arange(5,43,2)
signals_Fourier_1e6 = np.zeros((N_G1_F.size, index.size, number_of_bins))
n = 0
t_cpu_F_1e6 = np.zeros((N_G1_F.size, index.size))
for j in N_G1_F:
    for i in range(index.size):
        file_name = './Fourier_' + str(j) + '/TC_Fourier_NG1_' + str(j) + "_ind_" + str(i) + '_detector1'
        [width, height, simulation_data] = import_detector_signal_from_file(file_name)
        cropped_signal = simulation_data[1,(int(Xbins/2 - window_size_in_pixels/2)):(int(Xbins/2 + window_size_in_pixels/2))]
        rebinned_signal = np.mean(cropped_signal.reshape(number_of_bins,10),1)
        #plt.plot(rebinned_signal)
        #plt.show()
        signals_Fourier_1e6[n,i,:] = rebinned_signal / (np.sum(rebinned_signal))
        log_file = './Fourier_' + str(j) + '/TC_Fourier_NG1_' + str(j) + "_ind_" + str(i) + '.log'
        t_cpu_F_1e6[n,i] = get_simulation_time(log_file)[0]
    print(np.sum(signals_Fourier_1e6[n,:,:]))
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
    print(np.sum(signals_Fourier_4e6[n,:,:]))
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
    print(np.sum(signals_Fourier_4e6[n,:,:]))
    n += 1



#Wave propagator samples
all_signal_wp = np.load("/media/gerhard/SeagateExpansionDrive/PhD/MC/EGSnrc18/EGS_XGI_pub2/egs_xgi_home/Parameter_tests/TalbotCarpets/WP_TC_python/WP_rebinned_carpet.npy")
signal_WP = np.zeros((index.size, number_of_bins))
signal_WP = np.transpose(all_signal_wp[:,index])
for j in range(index.size):
    signal_WP[j,:] = signal_WP[j,:] / (np.sum(signal_WP[j,:]))
    #plt.plot(signal_WP[j,:])
    #plt.show()
    print((np.sum(signal_WP[j,:])))
print(np.sum(signal_WP))

#calculate difference measurements


#visulaize

mean_WP = np.mean(signal_WP)
mean_H = np.mean(signals_Huygens, (1,2))
mean_F_4e6 = np.mean(signals_Fourier_4e6, (1,2))

var_WP = np.var(signal_WP)
var_H = np.var(signals_Huygens, (1,2))
var_F_1e6 = np.var(signals_Fourier_1e6, (1,2))
var_F_4e6 = np.var(signals_Fourier_4e6, (1,2))
var_F_8e6 = np.var(signals_Fourier_8e6, (1,2))

print("var_WP: " + str(var_WP))
print("var_H: " + str(var_H))
print("var_F_1e6: " + str(var_F_1e6))
print("var_F_4e6: " + str(var_F_4e6))
print("var_F_8e6: " + str(var_F_8e6))


#calculate Pearson correlations:
cor_H = np.zeros(N_G1_H.size)
for j in range(N_G1_H.size):
    cov = np.cov(signal_WP.flatten(),signals_Huygens[j,:,:].flatten())
    #print(cov)
    cor_H[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
print("cor_H: " + str(cor_H))

cor_H_5e5 = np.zeros(N_G1_H.size)
for j in range(N_G1_H.size):
    cov = np.cov(signal_WP.flatten(),signals_Huygens_5e5[j,:,:].flatten())
    #print(cov)
    cor_H_5e5[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
print("cor_H_5e5: " + str(cor_H))


cor_F_8e6 = np.zeros(N_G1_F.size)
for j in range(N_G1_F.size):
        cov = np.cov(signal_WP.flatten(), signals_Fourier_8e6[j,:,:].flatten())
        cor_F_8e6[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
print("cor_F_8e6: " + str(cor_F_8e6))


cor_F_4e6 = np.zeros(N_G1_F.size)
for j in range(N_G1_F.size):
        cov = np.cov(signal_WP.flatten(), signals_Fourier_4e6[j,:,:].flatten())
        cor_F_4e6[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
print("cor_F_4e6: " + str(cor_F_4e6))

cor_F_1e6 = np.zeros(N_G1_F.size)
for j in range(N_G1_F.size):
        cov = np.cov(signal_WP.flatten(), signals_Fourier_1e6[j,:,:].flatten())
        cor_F_1e6[j] = cov[1,0] / np.sqrt(cov[0,0] * cov[1,1])
print("cor_F_1e6: " + str(cor_F_1e6))


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
print("get relative cpu time closest case: ")
d_corr_1e6 = (np.abs(cor_H[-1] / cor_F_1e6[:] - 1))
d_corr_1e6_H2 = (np.abs(cor_H[-2] / cor_F_1e6[:] - 1))
d_corr_4e6 = (np.abs(cor_H[-1] / cor_F_4e6[:] - 1))
d_corr_8e6 = (np.abs(cor_H[-1] / cor_F_8e6[:] - 1))

argmin_1e6 = np.argmin(d_corr_1e6)
argmin_1e6_H2 = np.argmin(d_corr_1e6_H2)
argmin_4e6 = np.argmin(d_corr_4e6)
argmin_8e6 = np.argmin(d_corr_8e6)


print("F_1e6:")
print(argmin_1e6)
print(d_corr_1e6[argmin_1e6])
print(cor_F_1e6[argmin_1e6])
print(N_G1_F[argmin_1e6])
print(mean_cpu_time_H[-1] / mean_cpu_time_F_1e6[argmin_1e6])

print("F_1e6_H2:")
print(argmin_1e6_H2)
print(d_corr_1e6_H2[argmin_1e6_H2])
print(cor_F_1e6[argmin_1e6_H2])
print(N_G1_F[argmin_1e6_H2])
print(mean_cpu_time_H[-2] / mean_cpu_time_F_1e6[argmin_1e6])


print("F_4e6:")
print(argmin_4e6)
print(d_corr_4e6[argmin_4e6])
print(cor_F_4e6[argmin_4e6])
print(N_G1_F[argmin_4e6])
print(mean_cpu_time_H[-1] / mean_cpu_time_F_4e6[argmin_4e6])

print("F_8e6:")
print(argmin_8e6)
print(d_corr_8e6[argmin_8e6])
print(cor_F_8e6[argmin_8e6])
print(N_G1_F[argmin_8e6])
print(mean_cpu_time_H[-1] / mean_cpu_time_F_8e6[argmin_8e6])

print("----------")
print("get relative cpu time closest case: ")
d_corr_1e6_5e5 = (np.abs(cor_H_5e5[-1] / cor_F_1e6[:] - 1))
d_corr_1e6_H2_5e5 = (np.abs(cor_H_5e5[-2] / cor_F_1e6[:] - 1))
d_corr_4e6_5e5 = (np.abs(cor_H_5e5[-1] / cor_F_4e6[:] - 1))
d_corr_8e6_5e5 = (np.abs(cor_H_5e5[-1] / cor_F_8e6[:] - 1))

argmin_1e6_5e5 = np.argmin(d_corr_1e6_5e5)
argmin_1e6_H2_5e5 = np.argmin(d_corr_1e6_H2_5e5)
argmin_4e6_5e5 = np.argmin(d_corr_4e6_5e5)
argmin_8e6_5e5 = np.argmin(d_corr_8e6_5e5)


print("F_1e6:")
print(argmin_1e6_5e5)
print(d_corr_1e6_5e5[argmin_1e6_5e5])
print(cor_F_1e6[argmin_1e6_5e5])
print(N_G1_F[argmin_1e6_5e5])
print(mean_cpu_time_H_5e5[-1] / mean_cpu_time_F_1e6[argmin_1e6_5e5])

print("F_1e6_H2:")
print(argmin_1e6_H2_5e5)
print(d_corr_1e6_H2_5e5[argmin_1e6_H2_5e5])
print(cor_F_1e6[argmin_1e6_H2_5e5])
print(N_G1_F[argmin_1e6_H2_5e5])
print(mean_cpu_time_H_5e5[-2] / mean_cpu_time_F_1e6[argmin_1e6_5e5])


print("F_4e6:")
print(argmin_4e6_5e5)
print(d_corr_4e6_5e5[argmin_4e6_5e5])
print(cor_F_4e6[argmin_4e6_5e5])
print(N_G1_F[argmin_4e6_5e5])
print(mean_cpu_time_H_5e5[-1] / mean_cpu_time_F_4e6[argmin_4e6_5e5])

print("F_8e6:")
print(argmin_8e6_5e5)
print(d_corr_8e6_5e5[argmin_8e6_5e5])
print(cor_F_8e6[argmin_8e6_5e5])
print(N_G1_F[argmin_8e6_5e5])
print(mean_cpu_time_H_5e5[-1] / mean_cpu_time_F_8e6[argmin_8e6_5e5])


print("closest case (restricted to the last two Huygens parameter choices)")
print("correlation coeff:")
print(cor_H[-2])
print(cor_F_4e6[6])
print(" relative correlation coeff")
print(cor_H[-2] / cor_F_4e6[6])
print(mean_cpu_time_H[-2] / mean_cpu_time_F_4e6[6])


print("mean cpu time")
print("H")
print(mean_cpu_time_H[-1])

print("F4e6")
print(mean_cpu_time_F_4e6[argmin_4e6])


print("F8e6")
print(mean_cpu_time_F_8e6[argmin_8e6])

plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
plt.plot(mean_cpu_time_H_5e5, cor_H_5e5, label=r'Huygens, $5\times10^5$', linewidth=2, color='k')
plt.plot(mean_cpu_time_H, cor_H, label=r'Huygens, $10^6$', linewidth=2, color='m')
plt.plot(mean_cpu_time_F_1e6, cor_F_1e6, label=r'Fourier, $10^6$', linewidth=2, color='green')
plt.plot(mean_cpu_time_F_4e6, cor_F_4e6, label=r'Fourier, $4\times10^6$', linewidth=2, color='tab:orange')
plt.plot(mean_cpu_time_F_8e6, cor_F_8e6, label=r'Fourier, $8\times10^6$', linewidth=2, color='blue')
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
plt.savefig("cor_tcpu_Huygens_Fourier_10carpetPositions_NH_ext_3.pdf", dpi=500, format="pdf")
plt.show()


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


#show comparison cases
plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
plt.plot([mean_cpu_time_F_1e6[0], mean_cpu_time_H[-1]],[cor_H[-1], cor_H[-1]], linewidth=2, color='grey', linestyle='dashed')
plt.plot([mean_cpu_time_F_1e6[0], mean_cpu_time_H[-1]],[cor_H_5e5[-1], cor_H_5e5[-1]],linewidth=2, color='grey', linestyle='dashed')
plt.plot(mean_cpu_time_H_5e5, cor_H_5e5, label="Huygens 5e5", linewidth=2, color='k')
plt.plot(mean_cpu_time_H, cor_H, label="Huygens 1e6", linewidth=2, color='m')
plt.plot(mean_cpu_time_F_1e6, cor_F_1e6, label="Fourier 1e6", linewidth=2, color='green')
plt.plot(mean_cpu_time_F_4e6, cor_F_4e6, label="Fourier 4e6", linewidth=2, color='tab:orange')
plt.plot(mean_cpu_time_F_8e6, cor_F_8e6, label="Fourier 8e6", linewidth=2, color='blue')
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
plt.savefig("cor_tcpu_Huygens_Fourier_10carpetPositions_NH_ext_4.pdf", dpi=500, format="pdf")
plt.show()



plt.figure(figsize=(11.0236, 8.26772), dpi = 109)
plt.plot(mean_cpu_time_H_5e5[1:], cor_H_5e5[1:], label="Huygens 5e5", linewidth=2, color='k')
plt.plot(mean_cpu_time_H, cor_H, label="Huygens 1e6", linewidth=2, color='m')
plt.plot(mean_cpu_time_F_1e6, cor_F_1e6, label="Fourier 1e6", linewidth=2, color='green')
plt.plot(mean_cpu_time_F_4e6, cor_F_4e6, label="Fourier 4e6", linewidth=2, color='tab:orange')
plt.plot(mean_cpu_time_F_8e6, cor_F_8e6, label="Fourier 8e6", linewidth=2, color='blue')
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
plt.savefig("cor_tcpu_Huygens_Fourier_10carpetPositions_NH_ext_rest.pdf", dpi=500, format="pdf")
plt.show()
