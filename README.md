# rDock

[![.github/workflows/build_matrix.yml](https://github.com/CBDD/rDock/actions/workflows/build_matrix.yml/badge.svg?branch=main)](https://github.com/CBDD/rDock/actions/workflows/build_matrix.yml)

- [Quick start guide](#quick-start-guide)
  - [Requirements](#requirements)
  - [Compilation](#compilation)
  - [Testing](#testing)
  - [Installation](#installation)
  - [Next steps](#next-steps)
- [rDock legacy version](#rdock-legacy-version)

## Quick start guide

This guide will help you compile and install rDock from source.  
If you're looking for precompiled binaries, you can find them in the [releases page](https://github.com/CBDD/rDock/releases), or go directly to the [latest release](https://github.com/CBDD/rDock/releases/latest).  
In order to install precompiled binaries, you can follow the same steps as the [installation section below](#installation), but you won't need to compile the project.

### Requirements

make sure the following requirements are installed and available:

- make
- a c++ compiler (g++ by default)
- popt and development headers (libpopt0 and libpopt-dev in ubuntu)
- git (optional if you download the code directly)

if you're running ubuntu, you can get all of them by running

```
sudo apt update && sudo apt install -y make git libpopt0 libpopt-dev g++
```

if you're running `macOS`, make sure to use gcc instead of clang. You can install gcc and popt using homebrew

```
brew install gcc popt
```

you can also check requirements for other officially supported distributions in the [Dockerfiles](https://github.com/CBDD/rDock/blob/main/.github/docker) used for CI

### Compilation

clone and compile the project

```
git clone https://github.com/CBDD/rDock
cd rDock
make
```

on `macOS`, you may need to run `export CXX=gcc` before running make in order to use gcc instead of clang. Here is an example for gcc-14 installed with homebrew:

```
export CXX=/opt/homebrew/Cellar/gcc/14.2.0_1/bin/g++-14
```

if you have multiple cores available, you may want to speed up the build process running make like this (replace 4 with the number of parallel processes you'd like to run)

```
make -j 4
```

for advanced compiling options, see the Makefile in this folder or run `make help`

### Testing

once the compilation process is finished, run tests to validate the built binaries

```
make test
```

### Installation

select the location for rDock binaries, library and development headers to be installed.

then set the PREFIX environment variable to point to this folder, for example ~/.local

```
PREFIX=~/.local make install
```

if PREFIX is not set, it will default to /usr, installing rDock for all users (you'll need sudo unless you're root):

```
sudo make install
```

make sure to add the installation folders to your PATH and LD_LIBRARY_PATH if necessary, and to set the RBT_ROOT environment variable:

```
export PATH=~/.local/bin:$PATH
export LD_LIBRARY_PATH=~/.local/lib:$LD_LIBRARY_PATH
export RBT_ROOT=~/.local/rDock
```

on `macOS`, LD_LIBRARY_PATH needs to be replaced by `DYLD_LIBRARY_PATH`

you may want to add these lines to your profile/configuration files like ~/.bashrc

### Next steps

If everything went well, you should be able to run rDock binaries like `rbcavity` or `rbdock`.  
You can find more information about the available tools in the [rDock documentation](https://rdock.github.io/documentation/), including several tutorials to get you started.

## rDock legacy version

with each new release, we try to improve the build system and the codebase, making it easier to maintain and to add new features, and the new features are not always fully compatible with the old ones.  
for this reason, the legacy version of rDock was frozen and will only receive critical bug fixes.  
you can find the [latest legacy release here](https://github.com/CBDD/rDock/releases/tag/v24.04.204-legacy)

