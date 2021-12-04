#!/bin/bash
cd ..
for angle in 1 5 10 30 60 100
do
cd /atlas/data19/zyni/cepc_dual1/condor_output1
hadd -f ${angle}.root electron_${angle}GeV_*.root
echo "${angle}.root"
done 