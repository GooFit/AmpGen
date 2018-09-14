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

All of the programs expect a file piped to stdin. They all support `--help` as well,
but unlike most command line programs, you must also run the program, and variable names and settings will be printed as
it runs (this is due to the fact that options are not pre-declared, but are declared where they are used).

### Generator

You can run a generator with:

```shell
./Generator MyOpts.opt --nEvents=1000 --output=output.root
```

This JIT compiles the model, and then generates `1000` events.

### ConvertToSourceCode

This produces the source code that you can compile (and also compiles it for you if you monitor the logs). This is used like this:

```shell
./ConvertToSourceCode MyOpts.opt --sourceFile=MyFile.cpp
```

You can compile this file (use the same commands given in the info statement printed out):

```shell
/usr/bin/c++ -Ofast -shared -rdynamic --std=c++14 -fPIC MyFile.cpp -o MyFile.so
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

