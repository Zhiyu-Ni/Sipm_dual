
# Sipm_dual

GEANT4-based simulation of a single bar for DualReadout purpose.

## Installing

Currently this only works on Umich cluster. Might need some modification if using on other clusters. 
```bash
source g4env.sh
cmake -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/10.1.p03/x86_64-slc6-gcc62-opt/lib64/Geant4-10.1.3
cmake --build .
```

## Single event running
For the interactive runing:

Run the events defined in run.mac
```bash
./CEPC_CaloTiming -c silicon_no_no.cfg -m run.mac -o test
# -c config file (files followed by .cfg are all config files)
# -m particle source and number
# -o output file
```

## Bulk running
use the generate_temp.sh in `./template`
```bash
cd template 
source generate_temp.sh
. submit_square.sh
```
