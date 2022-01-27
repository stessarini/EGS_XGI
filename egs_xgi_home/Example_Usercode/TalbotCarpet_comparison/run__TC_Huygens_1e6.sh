#!/bin/bash

###############################################################################
#
#   Bash script to run one set of Talbot carpet simulations
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

#if run individual scripts:
#export EGS_HOME=/path/to/egs_xgi_home/
#if not using default:
#export HEN_HOUSE=/path/to//HEN_HOUSE/
#export EGS_CONFIG=/path/to/config/file
#export my_machine=my_machine_name


n=1000000
i=0
i_max=10
j=3
j_max=7

cd ..


while [ $j -lt $j_max ]
do
  the_dir=./TalbotCarpet_comparison/Huygens_1e6_$j/
  i=0
  while [ $i -lt $i_max ]
  do
    file_name="TC_Hugens_1e6_exp_"$j"_ind_"$i
    output_file=$file_name".log"
    mv $the_dir"/"$file_name".egsinp" .
    rm $the_dir/$output_file
    echo $the_dir'/'$file_name
    ../bin/$my_machine/Example_Usercode -p 521icru -i $file_name -r RefractiveIndexFile -n $n > $the_dir/$output_file &
    if [[ $(( $i % 5 )) == 4 ]]
    then
      echo "wait"
      wait
    fi
    i=$[$i+1]
  done
  wait
  i=0
  echo "move:"
  while [ $i -lt $i_max ]
  do
    file_name="TC_Hugens_1e6_exp_"$j"_ind_"$i
    echo $the_dir'/'$file_name
    mv $file_name".egsinp" $the_dir
    i=$[$i+1]
  done
  j=$[$j+1]
  wait
done



exit 0
