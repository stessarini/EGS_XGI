#!/bin/bash

export EGS_HOME=/media/gerhard/SeagateExpansionDrive/PhD/MC/EGSnrc18/EGS_XGI_extension/egs_xgi_home/
#if not using default:
export HEN_HOUSE=/home/gerhard/MC/EGSnrc18/HEN_HOUSE/
export EGS_CONFIG=/home/gerhard/MC/EGSnrc18/HEN_HOUSE/specs/gcc750
export my_machine=gcc750

n=250000000

cd ..
the_dir=./PW_Cylinder_TestCases/
for i in {B..L..2}
do
  input=2mmx1mmFOV_PlaneWave_Case_$i
  echo $the_dir"/"$input".egsinp"
  mv $the_dir"/"$input".egsinp" .
  ../bin/$my_machine/Example_Usercode_Score_Energy -p 521icru_gi -i $input -r RefractiveIndexFile -n $n > $the_dir"/"$input".log"
  mv $input".egsinp" $the_dir
  wait
done
exit 0
