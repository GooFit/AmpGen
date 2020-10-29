#!/bin/bash

set -evx
# from https://stackoverflow.com/questions/55342122/conda-activate-on-travis-ci
export PATH="$DEPS_DIR/miniconda/bin:$PATH"
. $(conda info --root)/etc/profile.d/conda.sh
conda activate env_${TRAVIS_OS_NAME}

echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
echo "Building under OS: $TRAVIS_OS_NAME, CXX =$CXX"
echo "Directories: DEPS_DIR=$DEPS_DIR; TRAVIS_BUILD_DIR=$TRAVIS_BUILD_DIR"
mkdir -p build.conda
cd build.conda
cmake .. 
cmake --build . -- -j2

set +evx

cd ..

./build.conda/bin/AmpGen.exe options/example_b2kstarll.opt --CompilerWrapper::Verbose --nEvents 10000

