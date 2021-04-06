#!/bin/bash

. $(conda info --root)/etc/profile.d/conda.sh
conda activate env_${TRAVIS_OS_NAME}

echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
echo "Building under OS: $TRAVIS_OS_NAME"

mkdir -p build
cd build
echo "CMake-ing, CXX = $CXX"
cmake .. -DCMAKE_CXX_COMPILER=clang -DUSE_SIMD=0 -DUSE_OPENMP=0 -DUSE_MVEC=0 -DENABLE_INSTALL=0
echo "Building ..."
cmake --build . -- -j2 
cd ..
echo "Running test job ..."
./build/bin/AmpGen.exe options/example_b2kstarll.opt --CompilerWrapper::Verbose --nEvents 10000
