# source /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.3/x86_64-slc6/setup.sh
# source /cvmfs/geant4.cern.ch/geant4/10.5/x86_64-slc6-gcc63-opt-MT/CMake-setup.sh
# export CXX=`which g++`
# export CC=`which gcc`
# export PATH=$PATH:/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.11.1/Linux-x86_64/bin
# source /cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc62-opt/setup.sh 
# source $ROOTSYS/bin/thisroot.sh
# export LIBRARY_PATH=/home/eno/dualReadout/fakelib:$LIBRARY_PATH
# export LD_LIBRARY_PATH=/home/eno/dualReadout/fakelib:$LD_LIBRARY_PATH

# set up ROOT

export SMDT_DIR=`pwd`
export ALRB_rootVersion=6.10.06-x86_64-slc6-gcc62-opt
LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export ATLAS_LOCAL_ROOT_BASE=$LOCAL_ROOT_BASE
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalROOTSetup.sh --rootVersion ${ALRB_rootVersion}

# set up CMAKE

lsetup cmake
lsetup "gcc gcc620_x86_64_slc6"

# set up GEANT4

source /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.2.0/x86_64-slc6-gcc62-opt/setup.sh                  # set up compiler
source /cvmfs/geant4.cern.ch/geant4/10.1.p03/x86_64-slc6-gcc62-opt/CMake-setup.sh               # set up environment for Geant4
#source /cvmfs/sft.cern.ch/lcg/releases/XercesC/3.1.3-b3bf1/x86_64-slc6-gcc62-opt/XercesC-env.sh # set up GDML reader
export CXX=`which g++`                                                                          # tell CMake about compiler used 
export CC=`which gcc`
export G4INC='/cvmfs/geant4.cern.ch/geant4/10.1.p03/x86_64-slc6-gcc62-opt/include/Geant4'
export USE_VISUALISATION=1





