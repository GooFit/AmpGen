echo -en 'travis_fold:start:script.build\\r'
echo "Building..."
set -evx

mkdir -p build
cd build
cmake ..
cmake --build . -- -j2

set +evx
echo -e 'travis_fold:end:script.build\\r'
echo -en 'travis_fold:start:script.test\\r'
echo "Testing..."
set -evx

ctest --output-on-failure

