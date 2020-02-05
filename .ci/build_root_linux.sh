
# if [["$TRAVIS_OS_NAME" = "osx" ]] ; then return; fi

set -evx

pushd $DEPS_DIR 

# root_v6.12.06.Linux-ubuntu16-x86_64-gcc5.4.tar.gz
# ROOT_URL="https://root.cern.ch/download/root_v6.12.06.Linux-ubuntu16-x86_64-gcc5.4.tar.gz"
# ROOT_URL="https://root.cern.ch/download/root_v6.13.08.Linux-fedora28-x86_64-gcc8.0.tar.gz"
# 
# if [[ ! -f "${DEPS_DIR}/root/bin/root-config" ]] ; then
#   echo "Downloading Root"
#   mkdir -p root
#   travis_retry wget --no-check-certificate --quiet -O - "${ROOT_URL}" | tar --strip-components=1 -xz -C root
# fi
wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --add channels conda-forge
conda install --quiet --yes -c conda-forge/label/gcc8 root 

source "${DEPS_DIR}/root/bin/thisroot.sh"
popd

set +evx
