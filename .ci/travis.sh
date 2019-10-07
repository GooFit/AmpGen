echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
set -evx

mkdir -p build
cd build
cmake .. -DUSE_OPENMP=0
cmake --build . -- -j2

set +evx

cd ..

./build/bin/Generator options/example_b2kstarll.opt --CompilerWrapper::Verbose --nEvents 1000


# echo -e 'travis_fold:end:script.build\\r'
# echo -en 'travis_fold:start:script.test\\r'
# echo "Testing..."
# set -evx

# ctest --output-on-failure

