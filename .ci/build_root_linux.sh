pushd $DEPS_DIR 

if [[ $1 == "osx" ]] ; then 
  wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
elif [[ $1 == "linux" ]] ; then 
  wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
fi

bash miniconda.sh -b -p $DEPS_DIR/miniconda
export PATH="$DEPS_DIR/miniconda/bin:$PATH"
hash -r
conda config --add channels conda-forge
conda config --set channel_priority strict

# conda install --quiet --yes -c conda-forge/label/gcc8 root_base doxygen
conda create --quiet --yes -n my_root_env root doxygen zstd=1.3.7 -c conda-forge
# conda init bash
# source "$DEPS_DIR/miniconda/bin/thisroot.sh"
# export CXX="$DEPS_DIR/miniconda/bin/g++"
popd
