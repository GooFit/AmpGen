# AmpGen

## In the LHCb framework

AmpGen builds as part of Gauss in the LHCb framework.

## Local builds

AmpGen can also be built locally. The procedure is the standard one for CMake packages:

```shell
mkdir build
cd build
cmake ..
make
```

Once built, there will be several programs and option files in the `bin` directory.

### CentOS7

In order to build stand-alone on CentOS7, you will need a valid development environment; the following line will work:

```shell
lb-run ROOT $SHELL
```

### LLVM

You can also build AmpGen with LLVM. The only change you might want when using Apple LLVM
is to specifically specify the location of the build tool for AmpGen's JIT:

```shell
-DAMPGEN_CXX=$(which c++)
```

### Developers:
If using a full LLVM stack with access to LLVM tools, you can use them on
the AmpGen codebase. (For example, on Apple you'll need to run `brew install llvm`; the
built-in Apple LLVM does not have these tools.

#### Clang-tidy
To use clang-tidy to fix all warnings and style modernizations:

Add:

```shell
-DCMAKE_CXX_CLANG_TIDY="clang-tidy;-fix"
```
to the cmake configure line, and *do not* build in parallel! `make -j1`

#### Clang-format

You can make sure the code base follows the LHCb style requirements. Just run:

```shell
git ls-files '*.cpp' '*.h' | xargs clang-format -style=file -i
```

in the main AmpGen directory to process all the files.

## Include what you use

I would recommend a docker image:

```bash
docker run --rm -it tuxity/include-what-you-use:clang_4.0
```

And, you'll need to build root:

```bash
apt update && apt install libxpm-dev libxft-dev libxext-dev vim libgsl0-dev
git clone https://github.com/root-project/root.git root-src --branch=v6-12-06
mkdir root-build
cd root-buidl
cmake ../root-src -Dminuit2=ON -Dmathmore=ON
```

Then, you can run IWYU:

```bash
cmake .. -DCMAKE_CXX_INCLUDE_WHAT_YOU_USE="include-what-you-use"
make 2> iwyu.out
fix_includes.py < iwyu.out
```

You may need to convert `bits/shared_ptr.h` into `memory`.


## Using the stand-alone programs

All standalone programs can accept both options files and command line arguments. 
They also support `--help` to print help for key arguments to the program. 
This will also run the program, as arguments can be defined throughout each of the programs rather than all defined at the beginning. 

### Generator

The standalone generator for models can be used as

```shell
./Generator MyOpts.opt --nEvents=1000 --output=output.root
```

Which generates 1000 events of the model described in MyOpts.opt and saves them to output.root.

### ConvertToSourceCode

This produces source code to evalute the PDF, and normalises for use with other generators such as EvtGen, i.e. P(max) < 1. This can be used as

```shell
./ConvertToSourceCode MyOpts.opt --sourceFile=MyFile.cpp
```

This can then be a compiled to a shared library using

```shell
g++ -Ofast -shared -rdynamic --std=c++14 -fPIC MyFile.cpp -o MyFile.so
```

You can then check the file by opening it it Python3, using the utility class provided:

```python
from ampgen import FixedLib
lib = FixedLib('MyModel.so')
print(lib.matrix_elements[0]) # Print first matrix element

import pandas as pd
model = pd.read_csv('Input.csv', nrows=100_000)
fcn1 = lib.FCN_all(model)

# or, a bit slower, but just to show flexibility:
fcn2 = model.apply(lib.FCN, axis=1)
```

In this example, I've converted the output to csv to make it easy to import in Python without ROOT present. `.matrix_elements` gives you access to the parts of the PDF, with a `.amp` function that you can feed with P and E (`.vars` holds the current values).

### Fitter
### Debugger
### DataConverter

