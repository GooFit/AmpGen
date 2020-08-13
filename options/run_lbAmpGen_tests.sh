branch=mstahl_AmpGen

tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
echo "TMPDIR=$tmp_dir"
cd $tmp_dir

wget https://gitlab.cern.ch/lhcb/Gauss/-/archive/$branch/Gauss-${branch}.zip?path=Gen/LbAmpGen -O LbAmpGen.zip >> /dev/null
unzip LbAmpGen.zip >> /dev/null
top=Gauss-$branch-Gen-LbAmpGen

for model in $top/Gen/LbAmpGen/models/*.opt ; do 
  filename=$(basename $model)
  without_ext=${filename%.*}
  if [ $without_ext == "DtoKpipipi_v1" ] || 
     [ $without_ext == "DtopiKpipi_v1" ] || 
     [ $without_ext == "DtoKKpipi_v1"  ] ; then continue ; fi  # these are old models that kept for backwards compatability, dont expect to be able reproduce exactly
  mkdir -p build/$without_ext
  $AMPGENROOT/build/bin/ConvertToSourceCode $model --Output build/$without_ext/new.cpp  >> /dev/null
  g++ -Ofast -shared -rdynamic --std=c++11 -fPIC build/$without_ext/new.cpp -o build/$without_ext/new.so
  g++ -Ofast -shared -rdynamic --std=c++11 -fPIC $top/Gen/LbAmpGen/src/${without_ext}.cpp -o build/$without_ext/gaussUpdate.so
  $AMPGENROOT/build/bin/lib_diff    $model --Lib=build/$without_ext/new.so --RefLib=build/$without_ext/gaussUpdate.so
done
