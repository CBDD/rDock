[![.github/workflows/build_matrix.yml](https://github.com/CBDD/rDock/actions/workflows/build_matrix.yml/badge.svg?branch=main)](https://github.com/CBDD/rDock/actions/workflows/build_matrix.yml)
# rDock

## Quick start guide (new build system)

### Requirements

make sure the following requirements are installed and available:

* make
* a c++ compiler (g++ by default)
* popt and development headers (libpopt0 and libpopt-dev in ubuntu)
* git (optional if you download the code directly)

if you're running ubuntu, you can get all of them by running

```
sudo apt update && sudo apt install -y make git libpopt0 libpopt-dev g++
```

you can also check requirements for other distributions in the [Dockerfiles](https://github.com/CBDD/rDock/blob/main/.github/docker) used for CI

### Compilation

clone and compile the project

```
git clone https://github.com/CBDD/rDock
cd rDock
make
```

if you have multiple cores available, you may want to speed up the build process running make like this (replace 4 with the number of parallel processes you'd like to run)  
```
make -j 4
```

for advanced compiling options, see the Makefile in this folder

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

you may want to add these lines to your profile/configuration files like ~/.bashrc

## Old build system (deprecated):

this build system has been removed. if you _really_ need to use it, you can find it in [the frozen legacy branch here](https://github.com/CBDD/rDock/releases/tag/v24.04.204-legacy)
