#!/bin/bash

export EGS_HOME=/media/gerhard/SeagateExpansionDrive/PhD/MC/EGSnrc18/EGS_XGI_extension/egs_xgi_home/
#if not using default:
export HEN_HOUSE=/home/gerhard/MC/EGSnrc18/HEN_HOUSE/
export EGS_CONFIG=/home/gerhard/MC/EGSnrc18/HEN_HOUSE/specs/gcc750
export my_machine=gcc750

n=1000000
i=0
i_max=10
j=3
j_max=7

cd ..


while [ $j -lt $j_max ]
do
  the_dir=./TalbotCarpet_comparison/Huygens_$j/
  i=0
  while [ $i -lt $i_max ]
  do
    file_name="TC_Hugens_exp_"$j"_ind_"$i
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
  while [ $i -lt $i_max ]
  do
    file_name="TC_Hugens_exp_"$j"_ind_"$i
    echo $the_dir'/'$file_name
    mv $file_name".egsinp" $the_dir
    i=$[$i+1]
  done
  j=$[$j+1]
  wait
done



exit 0
