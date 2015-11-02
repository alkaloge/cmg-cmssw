#!/bin/bash
export SCRAM_ARCH=slc6_amd64_gcc481
#following suggestion on https://savannah.cern.ch/bugs/?100688 to access eos in multiple instances
export XRD_ENABLEFORKHANDLERS=1
WORK=$1; shift
SRC=$1; shift
cd $SRC; 
eval $(scramv1 runtime -sh);
cd $WORK;
exec $*
