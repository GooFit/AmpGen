pushd $DEPS_DIR 

wget -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $DEPS_DIR/miniconda
export PATH="$DEPS_DIR/miniconda/bin:$PATH"
hash -r
conda config --add channels conda-forge
conda install --quiet --yes -c conda-forge root
. $DEPS_DIR/miniconda/etc/profile.d/conda.sh
conda activate root
#         . "/home/tim/miniconda/etc/profile.d/conda.sh"
#     else
#         export PATH="/home/tim/miniconda/bin:$PATH"

# source "$DEPS_DIR/miniconda/bin/thisroot.sh"
popd

