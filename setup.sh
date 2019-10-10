AMPGEN=$(pwd)
export AMPGENROOT=$(pwd)/
export AMPGEN=$(pwd)

export PATH=$AMPGENROOT/install/bin:$AMPGENROOT/kspipi/bin:$AMPGENROOT/bin:$PATH

export PATH=$(pwd)/install/bin:$(pwd)/kspipi/bin:$PATH

#export PYTHONPATH=$(pwd)/install/bin:$(pwd)/kspipi/bin:$PYTHONPATH
if [ $(hostname | grep lxplus) ]
then
    source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_94python3 x86_64-centos7-gcc8-opt
fi
