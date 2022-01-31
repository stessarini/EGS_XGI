###############################################################################
#
#   Python analysis for PW_Cylinder image simulations
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


from import_detector_signal import get_detector_dimensions, get_number_of_histories, import_detector_signal_from_file, calculate_phase_stepping_curve



#number of pixels of the final images
nPixelX = 100
nPixelY = 50

#medium constants
mu_SI = 9.8542424227788743
mu_AIR = 0.00083461802611759059
mu_POLY = 0.40132757723110096

delta_SI = 1.2053424978669706e-06
delta_POLY = 5.8040247585592379e-07


print("############################################")
print("Import data...")

binary_file_name = '2mmx1mmFOV_PlaneWave_Reference_G12d0p5_0_detector1'
#dimensions of the detector
[FOV_x, FOV_y] = get_detector_dimensions(binary_file_name)
print("[FOV_x, FOV_y]: ")
print([FOV_x, FOV_y])
N_hist = get_number_of_histories(binary_file_name)
print("number of histories: " + str(N_hist))

print("load: " + binary_file_name)
[width, height, simulation_data_ref] = import_detector_signal_from_file(binary_file_name)

binary_file_name = '2mmx1mmFOV_PlaneWave_Reference_G12d0p5_0_NOP_detector1.bin'
print("load: " + binary_file_name)
simulation_data_NOP_ref = np.reshape(np.fromfile(binary_file_name,dtype=np.uintc),(height,width))



binary_file_name = '2mmx1mmFOV_PlaneWave_Cylinder_G12d0p5_0_detector1'
print("load: " + binary_file_name)
[width, height, simulation_data_sam] = import_detector_signal_from_file(binary_file_name)


binary_file_name = '2mmx1mmFOV_PlaneWave_Cylinder_G12d0p5_0_NOP_detector1.bin'
print("load: " + binary_file_name)
simulation_data_NOP_sam = np.reshape(np.fromfile(binary_file_name,dtype=np.uintc),(height,width))



total_signal_ref = np.sum(simulation_data_ref)
relative_singal_ref = total_signal_ref / N_hist





print("############################################")
print("Number of paths renormalization...")
#prevent divisions by 0
nonzero_ref = np.greater(simulation_data_NOP_ref, 0)
MC_signal_ref_NOP = np.zeros(simulation_data_ref.shape)
MC_signal_ref_NOP[nonzero_ref] = simulation_data_ref[nonzero_ref] / (simulation_data_NOP_ref[nonzero_ref]**2)
normalization_const_ref = total_signal_ref / np.sum(MC_signal_ref_NOP)
#keep total signal constant
MC_signal_ref_NOP = MC_signal_ref_NOP * normalization_const_ref
del nonzero_ref, simulation_data_ref, simulation_data_NOP_ref


total_signal_sam = np.sum(simulation_data_sam)
nonzero_sam = np.greater(simulation_data_NOP_sam, 0)
MC_signal_sam_NOP = np.zeros(simulation_data_sam.shape)
MC_signal_sam_NOP[nonzero_sam] = simulation_data_sam[nonzero_sam] / (simulation_data_NOP_sam[nonzero_sam]**2)
normalization_const_sam = total_signal_sam / np.sum(MC_signal_sam_NOP)
MC_signal_sam_NOP = MC_signal_sam_NOP * normalization_const_sam
del nonzero_sam, simulation_data_sam, simulation_data_NOP_sam





print("############################################")
print("Add primary and secondary signals...")

binary_file_name = '2mmx1mmFOV_PlaneWave_Reference_G12d0p5_0_SecondarySignal_detector1'
print("load: " + binary_file_name)
[width, height, simulation_data_sec_ref] = import_detector_signal_from_file(binary_file_name)
MC_signal_ref_NOP = MC_signal_ref_NOP + simulation_data_sec_ref
del simulation_data_sec_ref

binary_file_name = '2mmx1mmFOV_PlaneWave_Cylinder_G12d0p5_0_SecondarySignal_detector1'
print("load: " + binary_file_name)
[width, height, simulation_data_sec_sam] = import_detector_signal_from_file(binary_file_name)
MC_signal_sam_NOP = MC_signal_sam_NOP + simulation_data_sec_sam
del simulation_data_sec_sam




print("############################################")
print("Calculate phase stepping curves...")


#number of phase steps to perform
nPhaseSteps = 10
period_G2 = 0.0001
duty_cycle_G2 = 0.5 #other duty cycles currently not supported

#it is assumed that
#	the FOV is a integer multiple of period_G2/2
#	period_G2/2 is an integer multiple of data_point width

#nMask number of MC bins corresponding to the width of an absorbing section of G2
#i.e. how many consecutive bins are blocked by one slab of the grating
nMask = int(np.round(duty_cycle_G2 * width * period_G2 / FOV_x))
print("bins in x: " + str(width))
print("bins in y: " + str(height))
print("bins per pixel x: " + str(width / nPixelX))
print("bins per pixel y: " + str(height / nPixelY))
print("nMask " + str(nMask))

simulated_psc_ref = calculate_phase_stepping_curve(nPixelX, nPixelY, nMask, nPhaseSteps, MC_signal_ref_NOP)
simulated_psc_sam = calculate_phase_stepping_curve(nPixelX, nPixelY, nMask, nPhaseSteps, MC_signal_sam_NOP)

del MC_signal_ref_NOP, MC_signal_sam_NOP





print("############################################")
print("Get absorption and differenetial phase projections from phase stepping curves...")

absorption_signal = np.zeros([nPixelY, nPixelX])
differtial_phase_signal = np.zeros([nPixelY, nPixelX])

for nPixelCounterX in range(nPixelX):
	for nPixelCounterY in range(nPixelY):
		fft_Ref = np.fft.fft(simulated_psc_ref[:, nPixelCounterY, nPixelCounterX]);
		fft_Sam = np.fft.fft(simulated_psc_sam[:, nPixelCounterY, nPixelCounterX]);
		absorption_signal[nPixelCounterY, nPixelCounterX] = 1.0 - np.real(fft_Sam[0]/fft_Ref[0]);
		phi = np.angle(fft_Sam[1])-np.angle(fft_Ref[1]);
		if phi < -np.pi:
			differtial_phase_signal[nPixelCounterY, nPixelCounterX] = phi + np.ceil(-phi/(2.0*np.pi)) * 2.0 * np.pi;
		elif phi > np.pi:
			differtial_phase_signal[nPixelCounterY, nPixelCounterX] = phi - np.ceil(phi/(2.0*np.pi)) * 2.0 * np.pi;
		else:
			differtial_phase_signal[nPixelCounterY, nPixelCounterX] = phi;

del simulated_psc_ref, simulated_psc_sam
absorption_signal = np.flipud(absorption_signal)
absorption_signal = np.fliplr(absorption_signal)

differtial_phase_signal = np.flipud(differtial_phase_signal)
differtial_phase_signal = np.fliplr(differtial_phase_signal)







print("############################################")
print("Calculate Analytical absorption and differenetial phase profiles...")
#analytical values (using EGSnrc medium constants (ICRU512))
#cylinder radius
Radius=0.087;#[cm]

#pixel coordinates (centers)
pixel_width_cm = FOV_x / nPixelX
x_pixel_centers= np.linspace(-FOV_x/2.0 + pixel_width_cm / 2.0, FOV_x/2.0 - pixel_width_cm / 2.0, nPixelX);

#MC bins coordinates (centers)
dp_width_in_cm = FOV_x/width
x_data_point = np.linspace(-FOV_x / 2.0 + dp_width_in_cm / 2.0, FOV_x / 2.0 - dp_width_in_cm / 2.0 ,width)

#get projected cylinder thickness in the center of the bins
in_front_of_cylinder = np.less(x_data_point**2, Radius**2)
projected_cylinder_thickness = np.zeros(width)
projected_cylinder_thickness[in_front_of_cylinder] = 2.0 * np.sqrt(Radius**2 - x_data_point[in_front_of_cylinder]**2)

#Beer-Lambert law for given path lengths inside the medium
thy_abs_SI = np.exp(-mu_SI * projected_cylinder_thickness)
thy_abs_POLY = np.exp(-mu_POLY * projected_cylinder_thickness)
thy_abs_AIR = np.exp(-mu_AIR * projected_cylinder_thickness)

#get the absorption signal as an avereage over the bins within one pixel
reb_thy_abs_SI = np.mean(np.reshape(thy_abs_SI, (nPixelX, int(width/nPixelX))),axis=1)
reb_thy_abs_POLY = np.mean(np.reshape(thy_abs_POLY, (nPixelX, int(width/nPixelX))),axis=1)
reb_thy_abs_AIR = np.mean(np.reshape(thy_abs_AIR, (nPixelX, int(width/nPixelX))),axis=1)

#absorption profiles for SI or POLY cylinder
reb_thy_abs_signal_SI = 1.0 - reb_thy_abs_SI / reb_thy_abs_AIR
reb_thy_abs_signal_POLY = 1.0 - reb_thy_abs_POLY / reb_thy_abs_AIR



#Analytical differential phase signal
frac_talbot_dist=2.41966



thy_diff_phase_signal_SI = np.zeros(width);
thy_diff_phase_signal_POLY = np.zeros(width);

thy_diff_phase_signal_SI[in_front_of_cylinder] = (delta_SI * 2.0 * np.pi * frac_talbot_dist * 4.0 / period_G2) * x_data_point[in_front_of_cylinder] / (projected_cylinder_thickness[in_front_of_cylinder])
thy_diff_phase_signal_POLY[in_front_of_cylinder] = (delta_POLY * 2.0 * np.pi * frac_talbot_dist * 4.0 / period_G2) * x_data_point[in_front_of_cylinder] / (projected_cylinder_thickness[in_front_of_cylinder])

#remove infinities and NANs at edges

thy_diff_phase_signal_SI[np.isposinf(thy_diff_phase_signal_SI)]=np.pi
thy_diff_phase_signal_POLY[np.isposinf(thy_diff_phase_signal_POLY)]=np.pi

thy_diff_phase_signal_SI[np.isneginf(thy_diff_phase_signal_SI)]=-np.pi
thy_diff_phase_signal_POLY[np.isneginf(thy_diff_phase_signal_POLY)]=-np.pi


thy_diff_phase_signal_SI_reb = np.zeros(nPixelX);
thy_diff_phase_signal_POLY_reb = np.zeros(nPixelX);
in_front_of_cylinder = np.less(x_pixel_centers**2, Radius**2)
projected_cylinder_thickness = np.zeros(nPixelX)
projected_cylinder_thickness[in_front_of_cylinder] = 2.0 * np.sqrt(Radius**2 - x_pixel_centers[in_front_of_cylinder]**2)

thy_diff_phase_signal_SI_reb[in_front_of_cylinder] = (delta_SI * 2.0 * np.pi * frac_talbot_dist * 4.0 / period_G2) * x_pixel_centers[in_front_of_cylinder] / (projected_cylinder_thickness[in_front_of_cylinder])
thy_diff_phase_signal_POLY_reb[in_front_of_cylinder] = (delta_POLY * 2.0 * np.pi * frac_talbot_dist * 4.0 / period_G2) * x_pixel_centers[in_front_of_cylinder] / (projected_cylinder_thickness[in_front_of_cylinder])


thy_diff_phase_signal_SI_reb[np.greater(thy_diff_phase_signal_SI_reb, np.pi)] = np.pi
thy_diff_phase_signal_SI_reb[np.greater(-np.pi, thy_diff_phase_signal_SI_reb)] = -np.pi
thy_diff_phase_signal_POLY_reb[np.greater(thy_diff_phase_signal_POLY_reb, np.pi)] = -np.pi
thy_diff_phase_signal_POLY_reb[np.greater(-np.pi, thy_diff_phase_signal_POLY_reb)] = np.pi


#remove infinities and NANs at edges
thy_diff_phase_signal_SI[np.isinf(thy_diff_phase_signal_SI)]=0
thy_diff_phase_signal_POLY[np.isinf(thy_diff_phase_signal_POLY)]=0





print("############################################")
print("Set plot parameters...")
#
fig_width = 11.0236
fig_height = 8.26772
my_dpi = 109
#redefine in mm:
FOV_y = 1.0 #[mm]
#prepare lines to indicate the profile plots
x_coordinates = np.linspace(-1.0,1.0,10)
x_border_for_average_left = np.array([-1.0,-0.95])
x_border_for_average_right = np.array([0.95, 1.0])
y_coord_POLY = FOV_y/2.0 - FOV_y/nPixelY/2.0 - int(13) * FOV_y/nPixelY
y_coord_SI = FOV_y/2.0 - FOV_y/nPixelY/2.0 - int(36) * FOV_y/nPixelY
y_coordinates_POLY = np.ones(10) * y_coord_POLY
y_coordinates_SI = np.ones(10) * y_coord_SI

y_coordinate_range_for_average_abs_POLY_upper = np.ones(10) * (FOV_y/2.0 - FOV_y/nPixelY/2.0 - int(13 + 9) * FOV_y/nPixelY)
y_coordinate_range_for_average_abs_POLY_lower = np.ones(10) * (FOV_y/2.0 - FOV_y/nPixelY/2.0 - int(13 - 9) * FOV_y/nPixelY)
y_coordinate_range_for_average_abs_SI_upper = np.ones(10) * (FOV_y/2.0 - FOV_y/nPixelY/2.0 - int(36 + 9) * FOV_y/nPixelY)
y_coordinate_range_for_average_abs_SI_lower = np.ones(10) * (FOV_y/2.0 - FOV_y/nPixelY/2.0 - int(36 - 9) * FOV_y/nPixelY)



print("Absorption projection - indicate areas for average")
plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.imshow(np.clip(absorption_signal, a_min = 0.0, a_max = 1.0), cmap='gray',extent=[-1,1,-0.5,0.5], vmin=0, vmax=1.0, aspect="auto",interpolation='nearest')
#plt.plot(x_coordinates, y_coordinates_POLY, color='m', marker="", linestyle='solid',linewidth=2)
plt.plot(x_coordinates, y_coordinate_range_for_average_abs_POLY_upper, color='m', marker="", linestyle='dashed',linewidth=2)
plt.plot(x_coordinates, y_coordinate_range_for_average_abs_POLY_lower, color='m', marker="", linestyle='dashed',linewidth=2)
#plt.plot(x_coordinates, y_coordinates_SI, color='tab:orange', marker="", linestyle='solid',linewidth=2)
plt.plot(x_coordinates, y_coordinate_range_for_average_abs_SI_upper, color='tab:orange', marker="", linestyle='dashed',linewidth=2)
plt.plot(x_coordinates, y_coordinate_range_for_average_abs_SI_lower, color='tab:orange', marker="", linestyle='dashed',linewidth=2)
ax = plt.gca()
plt.text(0.0,0.9,'(a)', {'color': 'white', 'fontsize': 54},transform=ax.transAxes)
plt.xticks(np.around(np.linspace(-1,1,5),decimals=2), fontsize=50)
plt.yticks(fontsize=50)
plt.xlabel('x (mm)',fontsize=54)
plt.ylabel('y (mm)',fontsize=54)
#name='ABS_' + str(nPixelX) + 'x' + str(nPixelY) + '_nP'+ str(nPhaseSteps) + '_nM' + str(nMask) + '.png'
plt.subplots_adjust(bottom=0.19, top=0.96, right=0.94, left= 0.24)
#plt.savefig(name, bbox_inches = "tight")
plt.savefig("Absorption_projection_indicate_average.pdf", dpi=300, format="pdf")
plt.show()




#average the absorption signal over a few pixel rows
averaged_MC_absorption_signal_SI = np.mean(absorption_signal[(int(36) - 9):(int(2 * nPixelY/3) + 9),:], axis = 0)
averaged_MC_absorption_signal_POLY = np.mean(absorption_signal[(int(13) - 9):(int(nPixelY/3) + 9),:], axis = 0)


print('absorption profile SI and POLY')
plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.plot(10 * x_pixel_centers, 100 * reb_thy_abs_signal_POLY, color='k', marker="", linestyle='solid',linewidth=2)
plt.plot(10 * x_pixel_centers, 100 * reb_thy_abs_signal_SI, color='k', marker="", linestyle='dashed', linewidth=2)
plt.plot(10 * x_pixel_centers, 100 * averaged_MC_absorption_signal_POLY, color='m', marker=".", linestyle='none')
plt.plot(10 * x_pixel_centers, 100 * averaged_MC_absorption_signal_SI, color='tab:orange', marker=".", linestyle='none')
ax = plt.gca()
plt.text(0.0,0.9,'(b)', {'color': 'black', 'fontsize': 54},transform=ax.transAxes)
plt.xticks(np.around(np.linspace(-1,1,5),decimals=2),fontsize=50)
plt.yticks(np.arange(0,5) * 20,fontsize=50)
plt.grid()
plt.xlabel('x (mm)',fontsize=54)
plt.ylabel('Absorption (%)',fontsize=54)
#name='ABS_profile_SI_Poly_' + str(nPixelX) + 'x' + str(nPixelY) + '_nP'+ str(nPhaseSteps) + '_nM' + str(nMask) + '.png'
plt.subplots_adjust(bottom=0.19, top=0.96, right=0.96, left= 0.19)
#plt.savefig(name, bbox_inches = "tight")
plt.savefig("Averaged_absorption_profile.pdf", dpi=300, format="pdf")
plt.show()

print('absorption profile SI and POLY')
plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.plot(10 * x_pixel_centers, 100 * np.abs(reb_thy_abs_signal_POLY - averaged_MC_absorption_signal_POLY), color='m', marker="", linestyle='solid',linewidth=2, label='Poly')
plt.plot(10 * x_pixel_centers, 100 * np.abs(reb_thy_abs_signal_SI - averaged_MC_absorption_signal_SI), color='tab:orange', marker="", linestyle='dashed', linewidth=2, label='Si')
ax = plt.gca()
plt.text(0.0,0.9,'(c)', {'color': 'black', 'fontsize': 54},transform=ax.transAxes)
plt.xticks(np.around(np.linspace(-1,1,5),decimals=2),fontsize=50)
plt.yticks(np.arange(0,5) * 2,fontsize=50)
plt.grid()
plt.xlabel('x (mm)',fontsize=54)
plt.ylabel(r'$\left|\epsilon_{abs}\right| (\%)$',fontsize=54)
#plt.legend()
#name='ABS_profile_SI_Poly_' + str(nPixelX) + 'x' + str(nPixelY) + '_nP'+ str(nPhaseSteps) + '_nM' + str(nMask) + '.png'
plt.subplots_adjust(bottom=0.19, top=0.96, right=0.96, left= 0.23)
#plt.savefig(name, bbox_inches = "tight")
plt.savefig("abs_Error_Averaged_absorption_profile.pdf", dpi=300, format="pdf")
plt.show()




print("DPC projection - indicate areas for average")
plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.imshow(differtial_phase_signal, cmap='gray',extent=[-1,1,-0.5,0.5], vmin=-np.pi, vmax=np.pi, aspect="auto", interpolation='nearest')
#plt.plot(x_coordinates, y_coordinates_POLY, color='m', marker="", linestyle='solid',linewidth=2)
plt.plot(x_coordinates, y_coordinate_range_for_average_abs_POLY_upper, color='m', marker="", linestyle='dashed',linewidth=2)
plt.plot(x_coordinates, y_coordinate_range_for_average_abs_POLY_lower, color='m', marker="", linestyle='dashed',linewidth=2)
#plt.plot(x_coordinates, y_coordinates_SI, color='tab:orange', marker="", linestyle='solid',linewidth=2)
plt.plot(x_coordinates, y_coordinate_range_for_average_abs_SI_upper, color='tab:orange', marker="", linestyle='dashed',linewidth=2)
plt.plot(x_coordinates, y_coordinate_range_for_average_abs_SI_lower, color='tab:orange', marker="", linestyle='dashed',linewidth=2)
ax = plt.gca()
plt.text(0.0,0.9,'(d)', {'color': 'white', 'fontsize': 54},transform=ax.transAxes)
plt.xticks(np.around(np.linspace(-1,1,5),decimals=2), fontsize=50)
plt.yticks(fontsize=50)
plt.xlabel('x (mm)',fontsize=54)
plt.ylabel('y (mm)',fontsize=54)
#name='ABS_' + str(nPixelX) + 'x' + str(nPixelY) + '_nP'+ str(nPhaseSteps) + '_nM' + str(nMask) + '.png'
plt.subplots_adjust(bottom=0.19, top=0.96, right=0.94, left= 0.24)
#plt.savefig(name, bbox_inches = "tight")
plt.savefig("dpc_projection_indicate_average.pdf", dpi=300, format="pdf")
plt.show()




averaged_MC_dp_signal_SI = np.mean(differtial_phase_signal[(int(36) - 9):(int(36) + 9),:], axis = 0)
averaged_MC_dp_signal_POLY = np.mean(differtial_phase_signal[(int(13) - 9):(int(13) + 9),:], axis = 0)


print('diff phase signal SI and POLY')
plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.plot(10 * x_data_point, thy_diff_phase_signal_POLY,color='k', marker="", linestyle='solid',linewidth=2)
plt.plot(10 * x_data_point, thy_diff_phase_signal_SI,color='k', marker="", linestyle='dashed',linewidth=2)
plt.plot(10 * x_pixel_centers, averaged_MC_dp_signal_POLY,color='m', marker=".", linestyle='none')
plt.plot(10 * x_pixel_centers, averaged_MC_dp_signal_SI, color='tab:orange', marker=".", linestyle='none')
ax = plt.gca()
plt.text(0.0,0.9,'(e)', {'color': 'black', 'fontsize': 54},transform=ax.transAxes)
plt.xticks(np.around(np.linspace(-1,1,5),decimals=2),fontsize=50)#([2.74, 4.11,5.49,6.86,8.23])
plt.yticks(fontsize=50)
plt.grid()
plt.ylim((-3.5,3.5))
plt.xlabel('x (mm)',fontsize=54)
plt.ylabel('DPC (rad)',fontsize=54)
#name='DPC_profile_SI_Poly' + str(nPixelX) + 'x' + str(nPixelY) + '_nP'+ str(nPhaseSteps) + '_nM' + str(nMask) + '.png'
plt.subplots_adjust(bottom=0.19, top=0.96, right=0.96, left= 0.19)
plt.savefig("Averaged_differential_phase_profile.pdf", dpi=300, format="pdf")
plt.show()





print('diff phase signal SI and POLY')
plt.figure(figsize=(fig_width, fig_height), dpi = my_dpi)
plt.plot(10 * x_pixel_centers, np.abs(thy_diff_phase_signal_POLY_reb - averaged_MC_dp_signal_POLY),color='m', marker="", linestyle='solid',linewidth=2)
plt.plot(10 * x_pixel_centers, np.abs(thy_diff_phase_signal_SI_reb - averaged_MC_dp_signal_SI), color='tab:orange', marker="", linestyle='dashed',linewidth=2)

ax = plt.gca()
plt.text(0.0,0.9,'(f)', {'color': 'black', 'fontsize': 54},transform=ax.transAxes)
plt.xticks(np.around(np.linspace(-1,1,5),decimals=2),fontsize=50)#([2.74, 4.11,5.49,6.86,8.23])
plt.yticks(np.around(np.linspace(0,0.75,4),decimals=2),fontsize=50)
plt.grid()
#plt.ylim((0,0.75))
plt.xlabel('x (mm)',fontsize=54)
plt.ylabel(r'$\left|\epsilon_{DPC}\right|$ (rad)',fontsize=54)
#name='DPC_profile_SI_Poly' + str(nPixelX) + 'x' + str(nPixelY) + '_nP'+ str(nPhaseSteps) + '_nM' + str(nMask) + '.png'
plt.subplots_adjust(bottom=0.19, top=0.96, right=0.96, left= 0.23)
plt.savefig("abs_Error_Averaged_differential_phase_profile.pdf", dpi=300, format="pdf")
plt.show()







rmse_abs_SI = np.sqrt(np.mean((reb_thy_abs_signal_SI - averaged_MC_absorption_signal_SI)**2))
rmse_abs_POLY = np.sqrt(np.mean((reb_thy_abs_signal_POLY - averaged_MC_absorption_signal_POLY)**2))
rmse_dpc_SI = np.sqrt(np.mean((thy_diff_phase_signal_SI_reb - averaged_MC_dp_signal_SI)**2))
rmse_dpc_POLY = np.sqrt(np.mean((thy_diff_phase_signal_POLY_reb - averaged_MC_dp_signal_POLY)**2))

print("rmse_abs_SI: " + str(rmse_abs_SI))
print("rmse_abs_POLY: " + str(rmse_abs_POLY))
print("rmse_dpc_SI: " + str(rmse_dpc_SI))
print("rmse_dpc_POLY: " + str(rmse_dpc_POLY))
