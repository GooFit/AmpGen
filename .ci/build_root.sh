pushd $DEPS_DIR 

os=$1
if [[ $os == "osx" ]] ; then 
  wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda_${os}.sh
elif [[ $os == "linux" ]] ; then 
  wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda_${os}.sh
fi

bash miniconda_${os}.sh -b -p $DEPS_DIR/miniconda_${os}
export PATH="$DEPS_DIR/miniconda/bin:$PATH"
hash -r
conda config --add channels conda-forge
conda config --set channel_priority strict

# conda install --quiet --yes -c conda-forge/label/gcc8 root_base doxygen
conda create --quiet --yes -n env_${os} root doxygen zstd=1.3.7 -c conda-forge
# conda init bash
# source "$DEPS_DIR/miniconda/bin/thisroot.sh"
# export CXX="$DEPS_DIR/miniconda/bin/g++"
popd
