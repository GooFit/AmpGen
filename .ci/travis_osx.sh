echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
echo "Building under OS: $TRAVIS_OS_NAME"
mkdir -p build
cd build
echo "CMake-ing"
cmake .. #  -DUSE_OPENMP=0 -DCMAKE_CXX_COMPILER=clang
echo "Building ..."
cmake --build . -- -j2 
cd ..
echo "Running test job ..."
./build/bin/Generator options/example_b2kstarll.opt --CompilerWrapper::Verbose --nEvents 1000


# echo -e 'travis_fold:end:script.build\\r'
# echo -en 'travis_fold:start:script.test\\r'
# echo "Testing..."
# set -evx

# ctest --output-on-failure

