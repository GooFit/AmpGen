#!/bin/bash

# from https://stackoverflow.com/questions/55342122/conda-activate-on-travis-ci
export HOME=/home/tim/
export PATH="$HOME/miniconda/bin:$PATH"
. $(conda info --root)/etc/profile.d/conda.sh
conda activate env_${TRAVIS_OS_NAME}

echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
echo "Building under OS: $TRAVIS_OS_NAME, CXX =$CXX"
set -evx

mkdir -p build.conda
cd build.conda
cmake .. 
cmake --build . -- -j2

set +evx

cd ..

./build.conda/bin/Generator options/example_b2kstarll.opt --CompilerWrapper::Verbose --nEvents 10000

