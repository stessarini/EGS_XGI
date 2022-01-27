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
	#grating
	###################################
	:start geometry:
		name            = grating
		library         = egs_ndgeometry
		type            = EGS_XYZGeometry
		x-planes        = -0.005  -0.0048 -0.0046 -0.0044 -0.0042 -0.004  -0.0038 -0.0036 -0.0034 -0.0032 -0.003  -0.0028 -0.0026 -0.0024 -0.0022 -0.002  -0.0018 -0.0016 -0.0014 -0.0012 -0.001  -0.0008 -0.0006 -0.0004 -0.0002  0.0 0.0002  0.0004 0.0006  0.0008  0.001   0.0012  0.0014  0.0016  0.0018 0.002   0.0022  0.0024  0.0026  0.0028  0.003   0.0032 0.0034  0.0036  0.0038  0.004   0.0042  0.0044  0.0046 0.0048  0.005
		y-planes        = -0.000450 0.000450
		z-planes        = -0.0025715551 0.0
		:start media input:
			media = SI521ICRU VACUUM
			set medium = 1 1
			set medium = 3 1
			set medium = 5 1
			set medium = 7 1
			set medium = 9 1
			set medium = 11 1
			set medium = 13 1
			set medium = 15 1
			set medium = 17 1
			set medium = 19 1
			set medium = 21 1
			set medium = 23 1
			set medium = 25 1
			set medium = 27 1
			set medium = 29 1
			set medium = 31 1
			set medium = 33 1
			set medium = 35 1
			set medium = 37 1
			set medium = 39 1
			set medium = 41 1
			set medium = 43 1
			set medium = 45 1
			set medium = 47 1
			set medium = 49 1
		:stop media input:
	:stop geometry:
	############################################

	############################################
	# Envelope containing everything
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
		inscribed geometries = grating
	:stop geometry:
	############################################

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
		type = HuygensZPlane
		position = 0.0
		x-limits = -0.0051 0.0051
		y-limits = -0.00000001 0.00000001
		angle limits = -$phi $phi
		number of paths = $NG1
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
		save at = ./TalbotCarpet_comparison/Huygens_1e6_$exp_/
		allowed to overwrite = yes
	 :stop optical element:
:stop splitting algorithm definition:
''')

numberOfDetectors = 500
zCoordinate_values = np.linspace(3.22621786, 9.67865358, numberOfDetectors)
index = np.arange(0,10) * 49
#print(zCoordinate_values[index])
exp = np.arange(0,7)
angles = np.arctan(0.0102 / zCoordinate_values[index])
#print(angles)
for j in range(exp.size):
	for i in range(index.size):
		egs_input_file_name = './Huygens_1e6_' + str(j) + '/TC_Hugens_1e6_exp_' + str(j) + "_ind_" + str(i) + '.egsinp'
		N_G1=2**exp[j] * 225
		d = inputfile_template.substitute(zCoordinate=str(zCoordinate_values[index[i]]), Seed1='97', Seed2='42', NG1=str(N_G1), exp_=str(exp[j]), phi=str(angles[i]))
		egs_input_file = open(egs_input_file_name, "w")
		egs_input_file.write(d)
		egs_input_file.close()
