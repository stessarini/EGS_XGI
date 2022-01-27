# -*- coding: utf-8 -*-
###############################################################################
#
#   python script to generate one set of Talbot carpet simulation input files
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
from string import Template
import numpy as np

inputfile_template=Template('''\
###############################################
###############################################
##					     ##
##		Geometry		     ##
##					     ##
###############################################
###############################################

:start geometry definition:
	###################################
	###################################
		# Auxiliary grating wafer with medium
		:start geometry:
			library = egs_ndgeometry
			type = EGS_XYZGeometry
			name = Grating_dummy
			x-planes= -0.005150 0.005150
			y-planes= -0.000450 0.000450
			z-planes= -0.02720174 0.0
			:start media input:
				media = SI521ICRU
				set medium = 0 0
			:stop media input:
		:stop geometry:

		# Grating wafer used in the simulation
		:start geometry:
			library = egs_ndgeometry
			type = EGS_XYZGeometry
			name = VACWAFER
			x-planes= -0.005100 0.005100
			y-planes= -0.000400 0.000400
			z-planes= -0.02720174 0.0
			:start media input:
				media = VACUUM
				set medium = 0 0
			:stop media input:
		:stop geometry:

	############################################

	############################################
	# start Envelope containing everything
	###########################################

	:start geometry:
		library = egs_ndgeometry
		type = EGS_XYZGeometry
		name = base
		x-planes= -0.005200 0.005200
		y-planes= -0.000500 0.000500
		z-planes= -0.05720174 $zCoordinate
		:start media input:
			media = VACUUM
			set medium = 1 0
		:stop media input:
	:stop geometry:
	:start geometry:
		library = egs_genvelope
		name = simuGeom
		base geometry = base
		inscribed geometries = VACWAFER
	:stop geometry:
	###########################################
	# stop Envelope containing everything
	###########################################

	simulation geometry = simuGeom
:stop geometry definition:
###############################################
###############################################
##					     ##
##		source           	     ##
##					     ##
###############################################
###############################################
:start source definition:
	:start source:
		library = egs_parallel_beam
		name = PerfectPlaneWaveSource
		:start shape:
			library = egs_rectangle
			rectangle = -0.005100 -0.00000001 0.005100 0.00000001
		:stop shape:
		direction 0 0 1
		charge = 0
		:start spectrum:
			type = monoenergetic
			energy = 0.020
		:stop spectrum:
	:stop source:

	 :start source:
		library = egs_transformed_source
		name = TransformedPlaneWaveSource
		source name = PerfectPlaneWaveSource
		:start transformation:
			translation = 0 0 -0.05720174
		:stop transformation:
	:stop source:

	simulation source = TransformedPlaneWaveSource
:stop source definition:

###############################################
###############################################
##					     ##
##		Transport parameters         ##
##					     ##
###############################################
###############################################

:start MC transport parameter:

	 Global ECUT = 1.0
	# Global PCUT = 0.02

:stop MC transport parameter:


###############################################
###############################################
##					     ##
##		run control		     ##
##					     ##
###############################################
###############################################

:start run control:
	 ncase = 3000000#100000
:stop run control:
###############################################
###############################################
##					     ##
##		Random Number Generator      ##
##					     ##
###############################################
###############################################
:start rng definition:
	type = ranmar
	initial seeds = 97 42
:stop rng definition:

###############################################
###############################################
##					     ##
##		Splitting algorithm	     ##
##                      		     ##
###############################################
###############################################
:start splitting algorithm definition:
	:start splitting algorithm:
		type = DefaultSplittingAlgorithm
		setup = PhaseGrating detector1
		source coherence = coherent # incoherent coherent
	:stop splitting algorithm:

	:start optical element:
		name = PhaseGrating
		type = FourierSeriesZPlane
		position = 0.0
		phase shift at design energy = 1.0 #in units of pi
		duty cycle = 0.5
		periodicity = 0.0004 #in [cm]
		highest fourier order = $NG1
		design energy = 0.020
		Epsilon = 0.000001
		model attenuation = 0 #0 false; 1 true
		norm at design energy = 1.0 1.0 #if constant set for a , b sections: >= 0
		direction correction = 0 #0 false; 1 true
	:stop optical element:

	:start optical element:
		name = detector1
		type = DetectorZPlane
		total signal pixels = 30600 3
		position = -0.0051 -0.00000001 $zCoordinate
		height = 0.00000002
 		width = 0.0102
 		write to binary = 1
 		write to pgm = 0
		count number of paths = yes
		save at = ./TalbotCarpet_comparison/Fourier_1e6_$NG1/
		allowed to overwrite = yes
	 :stop optical element:
:stop splitting algorithm definition:
''')

numberOfDetectors = 500
zCoordinate_values = np.linspace(3.22621786, 9.67865358, numberOfDetectors)
index = np.arange(0,10) * 49
#print(zCoordinate_values[index])
N_G1=np.arange(5,25,2)
for j in N_G1:
	for i in range(index.size):
		egs_input_file_name = './Fourier_1e6_' + str(j) + '/TC_Fourier_1e6_NG1_' + str(j) + "_ind_" + str(i) + '.egsinp'
		d = inputfile_template.substitute(zCoordinate=str(zCoordinate_values[index[i]]), Seed1='97', Seed2='42', NG1=str(j))
		egs_input_file = open(egs_input_file_name, "w")
		egs_input_file.write(d)
		egs_input_file.close()
