#!/bin/bash

export EGS_HOME=/media/gerhard/SeagateExpansionDrive/PhD/MC/EGSnrc18/EGS_XGI_extension/egs_xgi_home/
#if not using default:
export HEN_HOUSE=/home/gerhard/MC/EGSnrc18/HEN_HOUSE/
export EGS_CONFIG=/home/gerhard/MC/EGSnrc18/HEN_HOUSE/specs/gcc750
export my_machine=gcc750

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
