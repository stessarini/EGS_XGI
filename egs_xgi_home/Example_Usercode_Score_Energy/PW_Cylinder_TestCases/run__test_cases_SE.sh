#!/bin/bash

###############################################################################
#
#   EGS_XGI Example_Usercode_Score_Energy run scripz
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
export my_machine=my_machine_name

n=250000000

cd ..
the_dir=./PW_Cylinder_TestCases/
for i in {1..3}
do
  input=2mmx1mmFOV_PlaneWave_Case_$i
  echo $the_dir"/"$input".egsinp"
  mv $the_dir"/"$input".egsinp" .
  ../bin/$my_machine/Example_Usercode_Score_Energy -p 521icru_gi -i $input -r RefractiveIndexFile -n $n > $the_dir"/"$input".log"
  mv $input".egsinp" $the_dir
  wait
done
exit 0
