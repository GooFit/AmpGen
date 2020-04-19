
wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $DEPS_DIR/miniconda
export PATH="$DEPS_DIR/miniconda/bin:$PATH"
hash -r
conda config --add channels conda-forge
conda config --set channel_priority strict
# conda install --quiet --yes -c conda-forge/label/gcc8 root_base doxygen
conda create --quiet --yes -n my_root_env root_base doxygen -c conda-forge/label/gcc8
conda activate my_root_env

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

