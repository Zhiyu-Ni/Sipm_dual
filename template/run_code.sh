#!/bin/bash
cd /atlas/data19/zyni/cepc_dual1
source /atlas/data19/zyni/cepc_dual1/g4env.sh
cd /atlas/data19/zyni/cepc_dual1

/atlas/data19/zyni/cepc_dual1/CEPC_CaloTiming -c $1 -m $2 -o $3
