pushd $DEPS_DIR 
set -evx

os=$1
if [[ $os == "osx" ]] ; then 
  wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda_${os}.sh
elif [[ $os == "linux" ]] ; then 
  wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda_${os}.sh
fi

bash miniconda_${os}.sh -b -p $DEPS_DIR/miniconda
export PATH="$DEPS_DIR/miniconda/bin:$PATH"
hash -r
conda config --add channels conda-forge
# conda config --set channel_priority strict

conda create --yes -n env_${os} root doxygen -c conda-forge

set +evx 
popd
