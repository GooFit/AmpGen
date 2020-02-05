echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
echo "Building under OS: $TRAVIS_OS_NAME, CXX =$CXX"
set -evx

mkdir -p build
cd build
cmake .. -DCMAKE_CXX_COMPILER=$CXX
cmake --build . -- -j2

set +evx

cd ..

./build/bin/Generator options/example_b2kstarll.opt --CompilerWrapper::Verbose --nEvents 1000


# echo -e 'travis_fold:end:script.build\\r'
# echo -en 'travis_fold:start:script.test\\r'
# echo "Testing..."
# set -evx

# ctest --output-on-failure

