#!/bin/bash
export PATH=$PWD/install/bin:$PATH
export AMPGENROOT=$PWD
if [ $(hostname | grep lxplus) ]
then
     source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_98python3 x86_64-centos7-gcc10-opt
##    source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_100_LHCB_6 x86_64-centos7-gcc10-opt
fi

