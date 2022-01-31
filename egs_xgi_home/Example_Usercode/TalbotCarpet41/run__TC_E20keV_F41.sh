#!/bin/bash

###############################################################################
#
#   Example run script for EGS_XGI applications
#		Runs all 500 Talbot carpet simulations.
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

export EGS_HOME=/path/to/egs_xgi_home/
#if EGSnrc environment variables are not set globally:
export HEN_HOUSE=/path/to//HEN_HOUSE/
export EGS_CONFIG=/path/to/config/file
export my_machine=set_my_machine_name

n=8000000
j_max=501

cd ..
the_dir=./TalbotCarpet41/
j=1
while [ $j -lt $j_max ]
do
  file_name="Carpet_G14d0p5_F41_20keV_Sim_"$j
  output_file="output_"$file_name".txt"
  mv $the_dir"/"$file_name".egsinp" .
  rm $output_file
  echo $file_name
  ../bin/$my_machine/Example_Usercode -p 521icru -i $file_name -r PolystyreneSpheres -n $n > $the_dir"/out_"$file_name".txt" &
  if [[ $(( $j % 5 )) == 0 ]]
  then
    wait
  fi
  j=$[$j+1]
done

#clean up
j=1
while [ $j -lt $j_max ]
do
  file_name="Carpet_G14d0p5_F41_20keV_Sim_"$j
  mv $file_name".egsinp" $the_dir
  j=$[$j+1]
done


exit 0
