#!/bin/bash
cd ..
rm -rf config1
if [ ! -d "config1" ]; then
   mkdir config1
fi
rm -rf condor_output1
path=$(pwd)
cd ./template
Particle="electron" 
#Particle="pion"
#pa="e-" 
#pa="pi-"
#tag=""
temp_name="silicon_no_no.cfg"
ex_name="submit_square.sh"
echo "" > $ex_name
for Particle in electron
do
if [[ "$Particle" == "electron" ]]; then
  pa="e-"
else
  pa="pi-"
fi

for number in 10
do
#echo "the particle is $Particle ($pa), the number is $number /n"
for energy in 5  #50 100 150 250 #1 2 5 10 20 30 40 50 60 70 80 90 100
do
for tag in  $(seq 0 1 100)  #0 1 2 3 4 5 6 7 8 9 
do
  source temp_mac.sh $pa $energy $number > $path/config1/run_${energy}GeV_${Particle}_N${number}_${tag}.mac
  source temp_jdl.sh $path ${Particle} $energy $number $tag $temp_name > $path/config1/condor-jobs_${energy}GeV_${Particle}_N${number}_${tag}.jdl
  echo "condor_submit $path/config1/condor-jobs_${energy}GeV_${Particle}_N${number}_${tag}.jdl" >> $ex_name
done
done
done
done
cd ..
if [ ! -d "condor_output1" ]; then
   mkdir condor_output1
fi
cd ./template
