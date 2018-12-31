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


