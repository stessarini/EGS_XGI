#!/bin/bash

export EGS_HOME=/media/gerhard/SeagateExpansionDrive/PhD/MC/EGSnrc18/EGS_XGI_extension/egs_xgi_home/
#if not using default:
export HEN_HOUSE=/home/gerhard/MC/EGSnrc18/HEN_HOUSE/
export EGS_CONFIG=/home/gerhard/MC/EGSnrc18/HEN_HOUSE/specs/gcc750
export my_machine=gcc750

n=25000

cd ..
the_dir=./LabSetup/
j=0

file_name="LabSetup"$j
output_file="out_"$file_name".txt"
mv $the_dir"/"$file_name".egsinp" .
rm $output_file
echo $file_name
../bin/$my_machine/Example_Usercode -p 521icru -i $file_name -r PolystyreneSpheres -n $n #> $the_dir"/out_"$file_name".txt"

file_name="LabSetup"$j
mv $file_name".egsinp" $the_dir



exit 0
