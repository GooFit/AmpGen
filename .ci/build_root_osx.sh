echo "Building ROOT for OSX [DIR=${DEPS_DIR}]"

pushd $DEPS_DIR 

# root_v6.12.06.Linux-ubuntu16-x86_64-gcc5.4.tar.gz
ROOT_URL="https://root.cern.ch/download/root_v6.12.06.macosx64-10.12-clang90.tar.gz"

if [ ! -f "${DEPS_DIR}/root/bin/root-config" ] ; then
  echo "Downloading Root"
  mkdir -p root
  travis_retry wget --no-check-certificate --quiet -O - "${ROOT_URL}" | tar --strip-components=1 -xz -C root
fi

source "${DEPS_DIR}/root/bin/thisroot.sh"
popd
