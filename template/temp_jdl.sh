#!/bin/sh

#define parameters which are passed in.
path=$1
particle=$2
energy=$3
event=$4
tag=$5
cfg=$6
#define the template.
#on_exit_hold =   (ExitCode != 0)
#on_exit_remove =   (ExitCode == 0)
#job_machine_attrs = Machine
#job_machine_attrs_history_length = 3

#Output = ${path}/run_${energy}GeV_${particle}_N${event}_${tag}.out
#Error =  ${path}/run_${energy}GeV_${particle}_N${event}_${tag}.err
#Log =    ${path}/run_${energy}GeV_${particle}_N${event}_${tag}.log

#Output = ${path}/condor_output/sce_\$(cluster)_\$(process).stdout
#Error =  ${path}/condor_output/sce_\$(cluster)_\$(process).stderr
#Log =    ${path}/condor_output/sce_\$(cluster)_\$(process).condor
#Requirements = TARGET.Machine =?= "r510-0-1.privnet"

cat  << EOF
getenv=True
Executable = ${path}/template/run_code.sh
Output = ${path}/condor_output1/run_${energy}GeV_${particle}_N${event}_${tag}.out
Error =  ${path}/condor_output1/run_${energy}GeV_${particle}_N${event}_${tag}.err
Log =    ${path}/condor_output1/run_${energy}GeV_${particle}_N${event}_${tag}.log
Arguments = ${path}/$cfg ${path}/config1/run_${energy}GeV_${particle}_N${event}_${tag}.mac ${path}/condor_output1/${particle}_${energy}GeV_N${event}_${tag}
Queue 
EOF
